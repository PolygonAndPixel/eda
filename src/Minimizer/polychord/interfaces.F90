!> This allows for a simple C interface...


module interfaces_module
    use utils_module, only: dp
    use iso_c_binding, only: c_ptr
    implicit none

    interface run_polychord
        module procedure run_polychord_full, run_polychord_no_prior
    end interface run_polychord

contains


    subroutine run_polychord_full(loglikelihood, prior, settings_in,dumper,context)
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
#ifdef MPI
        use mpi_module,               only: initialise_mpi, finalise_mpi
        use mpi,                      only: MPI_COMM_WORLD
#endif
        implicit none

        interface
            function loglikelihood(theta,phi,context)
                import :: dp
                import :: c_ptr
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
                type(c_ptr), intent(in), value    :: context
            end function loglikelihood
        end interface
        interface
            function prior(cube,context) result(theta)
                import :: dp
                import :: c_ptr
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
                type(c_ptr), intent(in), value    :: context
            end function prior
        end interface
        interface
            subroutine dumper(log_Z,log_ZError,ndead,n_llh_calls, RTI,settings,context)
                use settings_module, only: program_settings
                use run_time_type_module
                import :: c_ptr
                real(dp), intent(in), value             :: log_Z
                real(dp), intent(in), value             :: log_ZError
                real(dp), intent(in), value             :: ndead
                real(dp), intent(in), value             :: n_llh_calls
                ! real(dp), intent(in), value             :: n_accepted
                type(run_time_info), intent(in)         :: RTI
                type(program_settings), intent(in)      :: settings
                type(c_ptr), intent(in), value          :: context
            end subroutine dumper
        end interface

        ! integer, optional :: context
        type(c_ptr), intent(in), value    :: context

        type(program_settings),intent(in)    :: settings_in
        type(program_settings)               :: settings

        real(dp), dimension(4) :: output_info
#ifdef MPI
        call initialise_mpi
#endif
        call initialise_random()
        settings = settings_in
        call initialise_settings(settings)
#ifdef MPI
        output_info = NestedSampling(loglikelihood,prior,settings,MPI_COMM_WORLD,dumper,context)
        call finalise_mpi
#else
        output_info = NestedSampling(loglikelihood,prior,settings,0,dumper,context)
#endif

    end subroutine run_polychord_full


    !===================== INTERFACE ===============================================
    subroutine run_polychord_no_prior(loglikelihood, settings, dumper, context)
        use settings_module,          only: program_settings
        use loglikelihood_c,          only: no_prior
        implicit none
        interface
            function loglikelihood(theta,phi,context)
                import :: dp
                import :: c_ptr
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
                ! integer, optional :: context
                type(c_ptr), intent(in), value    :: context
            end function loglikelihood
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
                ! real(dp), intent(in), value             :: n_accepted
                type(run_time_info), intent(in)         :: RTI
                type(program_settings), intent(in)      :: settings
                type(c_ptr), intent(in), value          :: context
            end subroutine dumper
        end interface
        ! integer, optional :: context
        type(c_ptr), intent(in), value    :: context
        type(program_settings),intent(in)    :: settings  ! The program settings
        call run_polychord(loglikelihood,no_prior,settings,dumper,context)
    ! contains
        ! function prior(cube) result(theta)
        !     implicit none
        !     real(dp), intent(in), dimension(:) :: cube
        !     real(dp), dimension(size(cube))    :: theta
        !     theta = cube
        ! end function prior
    end subroutine run_polychord_no_prior










    subroutine polychord_c_interface( &
            nlive, num_repeats, nprior, do_clustering, &
            feedback, precision_criterion, max_ndead, &
            boost_posterior, posteriors, equals, &
            cluster_posteriors, write_resume, &
            write_paramnames, read_resume, write_stats, &
            write_live, write_dead, write_prior, &
            update_files, nDims, nDerived, base_dir, &
            file_root, nGrade, grade_frac, grade_dims, &
            c_loglikelihood_ptr, c_prior_ptr, c_dumper_ptr, context) &
            bind(c,name='polychord_c_interface')

        ! use iso_c_binding
        use loglikelihood_c
        use utils_module,             only: STR_LENGTH, convert_c_string
        use ini_module,               only: default_params
        use params_module,            only: param_type
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
        use read_write_module,        only: write_paramnames_file

        ! ~~~~~~~ Local Variable Declaration ~~~~~~~
        implicit none

        type(c_funptr), intent(in), value   :: c_loglikelihood_ptr
        type(c_funptr), intent(in), value   :: c_prior_ptr
        type(c_funptr), intent(in), value   :: c_dumper_ptr
        integer(c_int), intent(in), value   :: nlive
        integer(c_int), intent(in), value   :: num_repeats
        integer(c_int), intent(in), value   :: nprior
        logical(c_bool), intent(in), value  :: do_clustering
        integer(c_int), intent(in), value   :: feedback
        real(c_double), intent(in), value   :: precision_criterion
        integer(c_int), intent(in), value   :: max_ndead
        real(c_double), intent(in), value   :: boost_posterior
        logical(c_bool), intent(in), value  :: posteriors
        logical(c_bool), intent(in), value  :: equals
        logical(c_bool), intent(in), value  :: cluster_posteriors
        logical(c_bool), intent(in), value  :: write_resume
        logical(c_bool), intent(in), value  :: write_paramnames
        logical(c_bool), intent(in), value  :: read_resume
        logical(c_bool), intent(in), value  :: write_stats
        logical(c_bool), intent(in), value  :: write_live
        logical(c_bool), intent(in), value  :: write_dead
        logical(c_bool), intent(in), value  :: write_prior
        integer(c_int), intent(in), value   :: update_files
        integer(c_int), intent(in), value   :: nDims
        integer(c_int), intent(in), value   :: nDerived
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: base_dir
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: file_root
        ! integer context ! Actually a pointer to a custom C++ class
        type(c_ptr), intent(in), value    :: context

        integer(c_int), intent(in), value             :: nGrade
        real(c_double), intent(in), dimension(nGrade) :: grade_frac
        integer(c_int), intent(in), dimension(nGrade) :: grade_dims

        type(program_settings)    :: settings  ! The program settings

        type(param_type),dimension(:),allocatable :: params         ! Parameter array
        type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

        integer :: i_grade

        settings%nlive               = nlive
        settings%num_repeats         = num_repeats
        settings%nprior              = nprior
        settings%do_clustering       = do_clustering
        settings%feedback            = feedback
        settings%precision_criterion = precision_criterion
        settings%max_ndead           = max_ndead
        settings%boost_posterior     = boost_posterior
        settings%posteriors          = posteriors
        settings%equals              = equals
        settings%cluster_posteriors  = cluster_posteriors
        settings%write_resume        = write_resume
        settings%write_paramnames    = write_paramnames
        settings%read_resume         = read_resume
        settings%write_stats         = write_stats
        settings%write_live          = write_live
        settings%write_dead          = write_dead
        settings%write_prior         = write_prior
        settings%update_files        = update_files
        settings%nDims               = nDims
        settings%nDerived            = nDerived
        settings%base_dir            = convert_c_string(base_dir)
        settings%file_root           = convert_c_string(file_root)
        allocate(settings%grade_frac(nGrade),settings%grade_dims(nGrade))
        settings%grade_frac = grade_frac
        settings%grade_dims = grade_dims

        if(settings%write_paramnames) then
            params = default_params(settings%nDims,'theta','\theta')
            derived_params = default_params(settings%nDerived,'phi','\phi')
            call write_paramnames_file(settings,params,derived_params)
        end if
        call c_f_procpointer(c_loglikelihood_ptr, f_loglikelihood_ptr)
        call c_f_procpointer(c_prior_ptr, f_prior_ptr)
        call c_f_procpointer(c_dumper_ptr, f_dumper_ptr)
        call run_polychord(loglikelihood,prior,settings,dumper,context)

    end subroutine polychord_c_interface

end module interfaces_module
