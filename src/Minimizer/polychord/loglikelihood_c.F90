!> Addition for the C interface for IceCube's outdated compiler...

module loglikelihood_c
    use utils_module, only: dp
    use iso_c_binding
    implicit none

    interface
        function c_loglikelihood(theta,nDims,phi,nDerived,context) bind(c)
            use iso_c_binding
            integer(c_int), intent(in), value :: nDims
            integer(c_int), intent(in), value :: nDerived
            real(c_double), intent(in),  dimension(nDims) :: theta
            real(c_double), intent(out),  dimension(nDerived) :: phi
            real(c_double) :: c_loglikelihood
            ! integer context ! A pointer to a C++ class
            type(c_ptr), intent(in), value    :: context
        end function c_loglikelihood
    end interface
    interface
        subroutine c_prior(cube,theta,nDims,context) bind(c)
            use iso_c_binding
            integer(c_int), intent(in), value :: nDims
            real(c_double), intent(in),  dimension(nDims) :: cube
            real(c_double), intent(out), dimension(nDims) :: theta
            ! integer context ! A pointer to a C++ class
            type(c_ptr), intent(in), value    :: context
        end subroutine c_prior
    end interface
    interface
        subroutine c_dumper(log_Z,log_ZError,ndead,n_llh_calls, n_accepted, live_params,&
            n_cluster, llh_best_fit, llh_worst_fit, nPar, context) bind(c)
            use iso_c_binding
            integer(c_int), intent(in), value :: nPar
            real(c_double), intent(in), value :: log_Z
            real(c_double), intent(in), value :: log_ZError
            real(c_double), intent(in), value :: ndead
            real(c_double), intent(in), value :: n_llh_calls
            real(c_double), intent(in), value :: n_accepted
            real(c_double), intent(in), dimension(nPar) :: live_params
            integer(c_int), intent(in), value :: n_cluster
            real(c_double), intent(in), value :: llh_best_fit
            real(c_double), intent(in), value :: llh_worst_fit

            type(c_ptr), intent(in), value    :: context
        end subroutine c_dumper
    end interface

    procedure(c_loglikelihood), pointer :: f_loglikelihood_ptr
    procedure(c_prior), pointer         :: f_prior_ptr
    procedure(c_dumper), pointer        :: f_dumper_ptr

contains
    ! Hack due to IceCube's old fortran compiler.
    function no_prior(cube,context) result(theta)
        implicit none
        real(dp), intent(in), dimension(:) :: cube
        real(dp), dimension(size(cube))    :: theta
        ! integer, optional                  :: context
        type(c_ptr), intent(in), value    :: context
        theta = cube
    end function no_prior

    function loglikelihood(theta,phi,context)
        implicit none
        real(dp), intent(in),  dimension(:) :: theta
        real(dp), intent(out),  dimension(:) :: phi
        real(dp) :: loglikelihood

        real (c_double),dimension(size(theta)) :: c_theta
        integer (c_int)                        :: c_nDims
        real (c_double),dimension(size(phi))   :: c_phi
        integer (c_int)                        :: c_nDerived
        real (c_double)                        :: c_loglike
        ! integer context
        type(c_ptr), intent(in), value    :: context

        c_nDims = size(theta)
        c_nDerived = size(phi)
        c_theta = theta
        c_loglike = f_loglikelihood_ptr(c_theta,c_nDims,c_phi,c_nDerived,context)
        phi = c_phi
        loglikelihood = c_loglike

    end function loglikelihood

    function prior(cube, context) result(theta)
        implicit none
        real(dp), intent(in), dimension(:) :: cube
        real(dp), dimension(size(cube))    :: theta

        integer (c_int)                       :: c_nDims
        real (c_double),dimension(size(cube)) :: c_cube
        real (c_double),dimension(size(cube)) :: c_theta
        ! integer context
        type(c_ptr), intent(in), value    :: context

        c_nDims = size(cube)
        c_cube = cube
        call f_prior_ptr(c_cube,c_theta,c_nDims,context)
        theta = c_theta

    end function prior

    subroutine dumper(log_Z,log_ZError,ndead,n_llh_calls, RTI,settings,context)
        use settings_module, only: program_settings
        use run_time_type_module
        implicit none
        real(dp), intent(in), value             :: log_Z
        real(dp), intent(in), value             :: log_ZError
        real(dp), intent(in), value             :: ndead
        real(dp), intent(in), value             :: n_llh_calls
        type(run_time_info), intent(in)         :: RTI
        type(program_settings), intent(in)      :: settings
        type(c_ptr), intent(in), value          :: context

        real(dp)                                           :: best_fit
        real(dp)                                           :: best_fit_2
        integer                                            :: idx_fit_cluster, idx_fit_live
        real(dp), dimension(settings%grade_dims(1))        :: params_best_fit
        real(dp), dimension(settings%grade_dims(1))        :: params_best_fit_2
        real(dp)                                           :: worst_fit
        integer                                            :: i_live, i_cluster

        real(c_double)                                    :: c_log_Z
        real(c_double)                                    :: c_log_ZError
        real(c_double)                                    :: c_ndead
        real(c_double)                                    :: c_n_llh_calls
        real(c_double)                                    :: c_n_accepted
        real(c_double), dimension(settings%grade_dims(1)) :: c_live_params
        integer(c_int)                                    :: c_n_cluster
        real(c_double)                                    :: c_llh_best_fit
        real(c_double)                                    :: c_llh_worst_fit
        integer(c_int)                                    :: c_nPar

        c_log_Z = log_Z
        c_log_ZError = log_ZError
        c_ndead = ndead
        c_n_llh_calls = n_llh_calls
        c_live_params = RTI%params_best_fit
        c_n_cluster = RTI%ncluster
        c_llh_best_fit = RTI%best_fit
        c_llh_worst_fit = worst_fit
        c_n_accepted = RTI%n_accepted
        c_nPar = settings%grade_dims(1)

        call f_dumper_ptr(c_log_Z, c_log_ZError, c_ndead, c_n_llh_calls, c_n_accepted, &
            c_live_params, c_n_cluster, c_llh_best_fit, c_llh_worst_fit, c_nPar, &
            context)
    end subroutine dumper

end module loglikelihood_c
