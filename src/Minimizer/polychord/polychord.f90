module polychord

    ! ~~~~~~~ Loaded Modules ~~~~~~~
    use utils_module,             only: dp
    use settings_module,          only: program_settings
    use interfaces_module,        only: run_polychord
    use loglikelihood_module,     only: loglikelihood,prior,setup_loglikelihood,dumper
    use iso_c_binding,            only: c_ptr

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none
    type(program_settings)                    :: settings  ! The program settings
    contain

    subroutine polyRun(nDims, nDerived, nlive, num_repeats, do_clustering, &
                       precision_criterion, base_dir, write_resume, &
                       read_resume, write_live, write_dead, write_stats, &
                       equals, posterior, clurster_posteriors, feedback &
                       update_files, boost_posterior, context)

    integer         ,intent(in) :: nDims
    integer         ,intent(in) :: nDerived
    integer         ,intent(in) :: nlive
    integer         ,intent(in) :: num_repeats
    logical         ,intent(in) :: do_clustering
    real(dp)        ,intent(in) :: precision_criterion
    character(len=*),intent(in) :: base_dir
    character(len=*),intent(in) :: file_root
    logical         ,intent(in) :: write_resume
    logical         ,intent(in) :: read_resume
    logical         ,intent(in) :: write_live
    logical         ,intent(in) :: write_dead
    logical         ,intent(in) :: write_stats
    logical         ,intent(in) :: equals
    logical         ,intent(in) :: posteriors
    logical         ,intent(in) :: cluster_posteriors
    integer         ,intent(in) :: feedback
    integer         ,intent(in) :: update_files
    real(dp)        ,intent(in) :: boost_posterior
    type(c_ptr), intent(in), value    :: context

    ! Initialise the system settings
    settings%nDims         = nDims ! 20
    settings%nDerived      = nDerived ! 1
    settings%nlive         = nlive ! 500
    settings%num_repeats   = num_repeats ! settings%nDims*5
    settings%do_clustering = do_clustering ! .false.

    settings%precision_criterion = precision_criterion ! 1d-3

    settings%base_dir      = base_dir ! 'chains'
    settings%file_root     = file_root ! 'demo_gaussian'

    settings%write_resume  = write_resume ! .false.
    settings%read_resume   = read_resume ! .false.
    settings%write_live    = write_live ! .true.
    settings%write_dead    = write_dead ! .false.
    settings%write_stats   = write_stats ! .false.

    settings%equals        = equals ! .false.
    settings%posteriors    = posteriors ! .false.
    settings%cluster_posteriors = cluster_posteriors ! .false.

    settings%feedback      = feedback ! 1
    settings%update_files  = update_files ! settings%nlive

    settings%boost_posterior= boost_posterior ! 5.0_dp

    call setup_loglikelihood(settings)
    call run_polychord(loglikelihood, prior, settings, dumper, context)


end program PolyChord
