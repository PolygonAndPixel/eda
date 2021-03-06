module nested_sampling_module
    use utils_module, only: dp
    use iso_c_binding, only: c_ptr

#ifdef MPI
    use mpi_module, only: get_mpi_information,mpi_bundle,is_root,linear_mode,catch_babies,throw_babies,throw_seed,catch_seed,broadcast_integers,mpi_synchronise
#else
    use mpi_module, only: get_mpi_information,mpi_bundle,is_root,linear_mode
#endif

    implicit none


    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood,prior,settings,mpi_communicator,dumper,context) result(output_info)
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp,calc_similarity_matrix,swap_integers,cyc,time
        use read_write_module
        use feedback_module
        use run_time_module,   only: run_time_info,replace_point,calculate_logZ_estimate,calculate_covmats,delete_cluster,update_posteriors,delete_outermost_point
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: do_clustering
        use generate_module,   only: GenerateSeed,GenerateLivePoints,GenerateLivePointsFromSeed
#ifdef MPI
        use utils_module, only: logzero,normal_fb,stdout_unit
#endif

        implicit none

        ! Program Inputs
        ! ==============

        interface
            function loglikelihood(theta,phi,context)
                import :: dp
                import :: c_ptr
                real(dp), intent(in), dimension(:)  :: theta
                real(dp), intent(out), dimension(:) :: phi
                real(dp) :: loglikelihood
                ! integer, optional :: context
                type(c_ptr), intent(in), value    :: context
            end function
        end interface
        interface
            function prior(cube,context) result(theta)
                import :: dp
                import :: c_ptr
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
                ! integer, optional :: context
                type(c_ptr), intent(in), value    :: context
            end function
        end interface
        interface
            subroutine dumper(log_Z,log_ZError,ndead,n_llh_calls,RTI,settings,context)
                use settings_module, only: program_settings
                use run_time_module, only: run_time_info
                import :: dp
                import :: c_ptr
                real(dp), intent(in), value             :: log_Z
                real(dp), intent(in), value             :: log_ZError
                real(dp), intent(in), value             :: ndead
                real(dp), intent(in), value             :: n_llh_calls
                type(run_time_info), intent(in)         :: RTI
                type(program_settings), intent(in)      :: settings
                type(c_ptr), intent(in), value          :: context
            end subroutine dumper
        end interface

        !> Program settings
        type(program_settings), intent(in) :: settings

        !> MPI handle
        integer, intent(in) :: mpi_communicator

        ! Program Outputs
        ! ===============
        ! 1) log(evidence)
        ! 2) error(log(evidence))
        ! 3) ndead
        ! 4) number of likelihood calls
        ! 5) log(evidence) + log(prior volume)
        real(dp), dimension(4) :: output_info

        ! Pointer to a C++ class
        ! integer, optional :: context
        type(c_ptr), intent(in), value    :: context

        ! Local variables
        ! ===============

        ! The run time info
        ! ----------------
        ! very important, see src/run_time_info.f90
        type(run_time_info) :: RTI

        ! The number of repeats within each parameter speed to do
        integer, dimension(size(settings%grade_dims)) :: num_repeats
        integer, dimension(size(settings%grade_dims)) :: nlike      ! Temporary storage for number of likelihood calls
        integer, dimension(size(settings%grade_dims)) :: nlikesum   ! rolling average of number of loglikelihood calls


        ! Temporary variables
        ! -------------------
        ! Seed point to generate babies from
        real(dp), dimension(settings%nTotal) :: seed_point
        ! Cholesky matrix to generate babies from
        real(dp), dimension(settings%nDims,settings%nDims)   :: cholesky
        ! Loglikelihood bound to generate babies from
        real(dp) :: logL
        ! New-born baby points, created by slice sampling routine
        real(dp), allocatable, dimension(:,:) :: baby_points

        integer :: cluster_id ! Cluster identifier


        ! Logical Switches
        ! ----------------
        logical :: need_more_samples
        logical :: temp_logical

        ! MPI process variable
        ! --------------------
        type(mpi_bundle) :: mpi_information



#ifdef MPI
        ! MPI specific variables
        ! ----------------------
        integer                            :: i_slave       ! Slave iterator
        integer                            :: slave_id      ! Slave identifier
        integer, dimension(:), allocatable :: slave_cluster ! The cluster the slave is currently working on
        real(dp) :: time0,time1,slice_time,wait_time

        ! Slave switch
        ! ------------
        ! This prevents slaves delivering points to incorrect clusters after clustering
        ! has reorganised the cluster indices
        integer ::  slave_epoch
        integer ::  master_epoch
#endif


        ! A) Initialisation
        !    ==============
        ! MPI initialisation
        mpi_information = get_mpi_information(mpi_communicator)

#ifdef MPI
        allocate(slave_cluster(mpi_information%nprocs-1)) ! Allocate the slave arrays
        slave_cluster = 1                          ! initialise with 1

        ! slave switch
        slave_epoch=0
        master_epoch=0
#endif

        ! Rolling loglikelihood calculation
        nlikesum=0


        !-------------------------------------------------------!
        if(is_root(mpi_information)) call write_opening_statement(settings) !
        !-------------------------------------------------------!



        ! Check if we actually want to resume
        if ( settings%read_resume .and. resume_file_exists(settings) ) then

            ! Read the resume file on root
            if(is_root(mpi_information)) then
                call read_resume_file(settings,RTI)
                ! -------------------------------------------- !
                call write_resuming(settings%feedback)
                ! -------------------------------------------- !
            end if


        else

            ! Delete any existing files if we're going to be producing our own new ones
            if(is_root(mpi_information).and.settings%write_resume) call delete_files(settings)

            ! Intialise the run by setting all of the relevant run time info, and generating live points
            if (settings%generate_from_seed) then
                call GenerateLivePointsFromSeed(loglikelihood,prior,settings,&
                    RTI,mpi_information,context)
            else
                call GenerateLivePoints(loglikelihood,prior,settings,RTI,&
                    mpi_information, context)
            end if
            if(is_root(mpi_information).and.settings%write_prior) then
                call write_prior_file(settings,RTI)
            end if

            do while(is_root(mpi_information) .and. RTI%nlive(1) > settings%nlive )
                call delete_outermost_point(settings,RTI)
            end do

            ! Write a resume file (as the generation of live points can be intensive)
            if(is_root(mpi_information).and.settings%write_resume) then
                call write_resume_file(settings,RTI)
                call rename_files(settings,RTI)
            end if


        end if

        if(is_root(mpi_information)) then
            num_repeats = RTI%num_repeats
            call write_num_repeats(num_repeats,settings%feedback)
        end if
#ifdef MPI
        call broadcast_integers(num_repeats,mpi_information)
#endif
        allocate(baby_points(settings%nTotal,sum(num_repeats)))



        ! B) Main loop body
        !    ==============

        if(is_root(mpi_information)) then

            ! -------------------------------------------- !
            call write_started_sampling(settings%feedback) !
            ! -------------------------------------------- !

            ! Definitely need more samples than this
            need_more_samples = RTI%ndead<settings%max_ndead .or. settings%max_ndead<0

            do while ( need_more_samples )


                ! 1) Generate a new live point
                !    -------------------------
                ! Generate a seed point --- update this
                seed_point = GenerateSeed(settings,RTI,cluster_id)

                ! Choose the cholesky decomposition for the cluster
                cholesky = RTI%cholesky(:,:,cluster_id)

                ! Get the loglikelihood contour we're generating from
                logL = RTI%logLp(cluster_id)


                if(linear_mode(mpi_information)) then
                    ! Linear Mode
                    ! -----------

                    ! Generate a new set of points within the likelihood bound of the late point
                    baby_points = SliceSampling(loglikelihood,prior,settings,&
                        logL,seed_point,cholesky,nlike,num_repeats,RTI,context)
#ifdef MPI
                else
                    ! Parallel mode
                    ! -------------

                    ! Recieve any new baby points from any slave currently sending
                    slave_id = catch_babies(baby_points,nlike,slave_epoch,mpi_information)

                    ! and throw seeding information back to slave (true => keep going)
                    call throw_seed(seed_point,cholesky,logL,mpi_information,slave_id,master_epoch,.true.)

                    ! set cluster_id to be the cluster identity of the babies just recieved
                    ! (saved in slave_cluster from the last send) and set slave_cluster to
                    ! be the bound just sent off.
                    call swap_integers(cluster_id,slave_cluster(slave_id))

#endif
                end if

                ! Add the likelihood calls to our counter
                RTI%nlike = RTI%nlike + nlike
                nlikesum  = nlikesum  + nlike


                ! See if this point is suitable to be added to the arrays
#ifdef MPI
                if( linear_mode(mpi_information) .or. master_epoch==slave_epoch ) then
#endif
                    if(replace_point(settings,RTI,baby_points,cluster_id) ) then

                        ! Check to see if we need more samples
                        need_more_samples = more_samples_needed(settings,RTI)

                        ! Update the posterior array every nlive iterations
                        if( cyc(RTI%ndead,settings%nlive) ) call update_posteriors(settings,RTI)

                        ! Update the resume files every settings%update_resume iterations,
                        ! or at the end of the run
                        if( cyc(RTI%ndead,settings%update_files) .and. settings%update_files > 0 ) then
                            if(settings%write_resume)                  call write_resume_file(settings,RTI)
                            if(settings%write_live)                    call write_phys_live_points(settings,RTI)
                            if(settings%write_dead)                    call write_dead_points(settings,RTI)
                            if(settings%write_stats)                   call write_stats_file(settings,RTI,nlikesum)
                            if(settings%equals.or.settings%posteriors) call write_posterior_file(settings,RTI)
                            call rename_files(settings,RTI)
                        end if


                        if(delete_cluster(settings,RTI)) then
#ifdef MPI
                            master_epoch = master_epoch+1
#endif
                        end if! Delete any clusters as necessary

                        ! update the clustering and covariance matrices every nlive iterations
                        if( cyc(RTI%ndead,settings%nlive) ) then
                            !--------------------------------------------!
                            call write_intermediate_results(settings,RTI,nlikesum)
                            nlikesum=0
                            !--------------------------------------------!
                            if(settings%do_clustering) then

                                ! If we want to cluster on sub dimensions, then do this first
                                if(allocated(settings%sub_clustering_dimensions)) then
                                    if( do_clustering(settings,RTI,settings%sub_clustering_dimensions) )  then
#ifdef MPI
                                        master_epoch = master_epoch+1
#endif
                                    end if
                                end if

                                if( do_clustering(settings,RTI) )  then
#ifdef MPI
                                    master_epoch = master_epoch+1
#endif
                                end if
                            end if
                            call calculate_covmats(settings,RTI)
                        end if
                    end if
#ifdef MPI
                end if
#endif

            end do ! End of main loop body

            ! Clean up the remaining live points
            if(settings%write_resume)                  call write_resume_file(settings,RTI)

            do while(RTI%ncluster > 0)
                call delete_outermost_point(settings,RTI)
                temp_logical = delete_cluster(settings,RTI)
            end do

            call update_posteriors(settings,RTI)
            if(settings%write_live)                    call write_phys_live_points(settings,RTI)
            if(settings%write_stats)                   call write_stats_file(settings,RTI,nlikesum)
            if(settings%equals.or.settings%posteriors) call write_posterior_file(settings,RTI)
            if(settings%write_dead)                    call write_dead_points(settings,RTI)
            call rename_files(settings,RTI)

            ! Create the output array
            ! (1) log evidence
            ! (2) variance in the log evidence
            ! (3) Number of dead points
            ! (4) Number of slow likelihood calls
            ! (5) log(evidence * prior volume)
            call calculate_logZ_estimate(RTI,output_info(1),output_info(2))
            output_info(3) = RTI%ndead
            output_info(4) = RTI%nlike(1)
            call dumper(output_info(1),output_info(2),output_info(3),&
                output_info(4),RTI,settings,context)
            ! ------------------------------------------------------------ !
            call write_final_results(output_info,settings%feedback)        !
            ! ------------------------------------------------------------ !

            ! call dumper(output_info(1),output_info(2),output_info(3),&
            !     output_info(4),RTI%live,RTI%nlive,RTI%ncluster,&
            !     RTI%logLp,settings%nDims,RTI%i,context)







            ! C) Clean up
            !    ========

#ifdef MPI
            ! MPI cleanup
            ! -----------
            ! Kill off the final slaves.
            ! If we're done, then clean up by receiving the last piece of
            ! data from each node (and throw it away) and then send a kill signal back to it
            do i_slave=mpi_information%nprocs-1,1,-1

                ! Recieve baby point from slave slave_id
                slave_id = catch_babies(baby_points,nlike,slave_epoch,mpi_information)

                ! Add the likelihood calls to our counter
                RTI%nlike = RTI%nlike + nlike

                ! Send kill signal to slave slave_id (note that we no longer care about seed_point, so we'll just use the last one
                call throw_seed(seed_point,cholesky,logL,mpi_information,slave_id,master_epoch,.false.)

            end do


        else !(myrank/=root)

            ! These are the slave tasks
            ! -------------------------
            !
            ! This is considerably simpler than that of the master.
            ! All slaves do is:
            ! 1) recieve a seed point,cholesky decomposition and loglikelihood
            !    contour from the master
            ! 2) using the above, generate a new set of baby points
            ! 3) send the baby points and nlike back to the master.


            ! On the first loop, send a nonsense set of baby_points
            ! to indicate that we're ready to start receiving

            baby_points = 0d0                     ! Avoid sending nonsense
            baby_points(settings%l0,:) = logzero  ! zero contour to ensure these are all thrown away
            nlike = 0                             ! no likelihood calls in this round
            call throw_babies(baby_points,nlike,slave_epoch,mpi_information)
            wait_time = 0
            slice_time = 0
            time1 = time()




            ! 1) Listen for a seed point being sent by the master
            !    Note that this also tests for a kill signal sent by the master
            do while(catch_seed(seed_point,cholesky,logL,slave_epoch,mpi_information))
                time0 = time()
                ! 2) Generate a new set of baby points
                baby_points = SliceSampling(loglikelihood,prior,settings,logL,
                    seed_point,cholesky,nlike,num_repeats,RTI,context)


                wait_time = wait_time + time0-time1
                time1 = time()
                slice_time = slice_time + time1-time0


                ! 3) Send the baby points back
                call throw_babies(baby_points,nlike,slave_epoch,mpi_information)

            end do

            if(slice_time<wait_time) then
                if(settings%feedback<=normal_fb) write(stdout_unit,'("Slave",I3,": Inefficient MPI parallisation, I spend more time waiting than slicing ", E17.8, ">", E17.8 )') mpi_information%rank, wait_time,slice_time
            else
                if(settings%feedback<=normal_fb) write(stdout_unit,'("Slave",I3,": efficient MPI parallisation; wait_time/slice_time= ", E17.8 )') mpi_information%rank, wait_time/slice_time
            end if

#endif
        end if !(myrank==root / myrank/=root)

#ifdef MPI
        call mpi_synchronise(mpi_information)
#endif




    end function NestedSampling












    !> This function checks whether we've done enough to stop
    function more_samples_needed(settings,RTI)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info,live_logZ
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: RTI

        logical :: more_samples_needed

        ! Set it to default to true
        more_samples_needed = .true.

        ! If we've put a maximum number of iterations on the algorithm, then
        ! we'll stop if we've reached that number.
        ! If we don't want a maximum number of iterations, then max_ndead should
        ! be set negative or 0
        if(settings%max_ndead>0 .and. RTI%ndead >= settings%max_ndead) then
            more_samples_needed = .false.
            return
        end if

        ! If the evidence in the live points is less than precision_criterion %
        ! of the total accumulated evidence, then stop.
        if( live_logZ(settings,RTI) < log(settings%precision_criterion) + RTI%logZ )  then
            more_samples_needed = .false.
            return
        end if


    end function more_samples_needed







end module nested_sampling_module
