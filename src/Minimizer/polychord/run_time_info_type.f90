module run_time_type_module
    use utils_module, only: dp
    implicit none

    !> The run time information.
    !!
    !! This is what needs to be saved in order to resume a run.
    !! Bundling these all into the same type enables easy passing of data
    !! from one fuction to another
    type run_time_info

        !> Number of dead points
        integer :: ndead

        !> The number currently evolving clusters
        integer :: ncluster
        !> The number of dead clusters
        integer :: ncluster_dead

        !> Total number of likelihood calls
        integer,allocatable,dimension(:) :: nlike

        !> The number of repeats within each parameter speed to do
        integer,allocatable, dimension(:)           :: num_repeats

        !> The number of live points in each cluster
        integer, allocatable, dimension(:) :: nlive
        !> The number of phantom points in each cluster
        integer, allocatable, dimension(:) :: nphantom
        !> The number of weighted posterior points in each cluster
        integer, allocatable, dimension(:) :: nposterior
        !> The number of equally weighted posterior points in each cluster
        integer, allocatable, dimension(:) :: nequals


        integer, allocatable, dimension(:) :: nposterior_dead
        integer                            :: nposterior_global
        integer, allocatable, dimension(:) :: nequals_dead
        integer                            :: nequals_global

        !> Live points
        real(dp), allocatable, dimension(:,:,:) :: live
        !> Phantom points
        real(dp), allocatable, dimension(:,:,:) :: phantom
        !> Posterior stack
        real(dp), allocatable, dimension(:,:,:) :: posterior_stack
        !> The number of posterior points in each cluster in the posterior stack
        integer, allocatable, dimension(:) :: nposterior_stack


        !> weighted posterior points
        real(dp), allocatable, dimension(:,:,:) :: posterior
        real(dp), allocatable, dimension(:,:,:) :: posterior_dead
        real(dp), allocatable, dimension(:,:)   :: posterior_global

        !> Equally weighted posterior points
        real(dp), allocatable, dimension(:,:,:) :: equals
        real(dp), allocatable, dimension(:,:,:) :: equals_dead
        real(dp), allocatable, dimension(:,:)   :: equals_global

        !> Pure nested sampling points
        real(dp), allocatable, dimension(:,:)   :: dead

        !> Covariance Matrices
        real(dp), allocatable, dimension(:,:,:) :: covmat
        !> Cholesky decompositions
        real(dp), allocatable, dimension(:,:,:) :: cholesky


        !> Global evidence estimate
        real(dp) :: logZ
        !> Global evidence^2 estimate
        real(dp) :: logZ2


        !> Local volume estimate
        real(dp), allocatable, dimension(:)   :: logXp
        !> global evidence volume cross correlation
        real(dp), allocatable, dimension(:)   :: logZXp
        !> Local evidence estimate
        real(dp), allocatable, dimension(:)   :: logZp
        real(dp), allocatable, dimension(:)   :: logZp_dead
        !> Local evidence^2 estimate
        real(dp), allocatable, dimension(:)   :: logZp2
        real(dp), allocatable, dimension(:)   :: logZp2_dead
        !> local evidence volume cross correlation
        real(dp), allocatable, dimension(:)   :: logZpXp
        !> local volume cross correlation
        real(dp), allocatable, dimension(:,:) :: logXpXq

        !> Minimum loglikelihoods
        real(dp), allocatable, dimension(:) :: logLp
        !> The minimum loglikelihood point within each cluster
        integer,allocatable, dimension(:)           :: i

        !> Maximum weight
        real(dp), allocatable, dimension(:) :: maxlogweight
        real(dp), allocatable, dimension(:) :: maxlogweight_dead
        real(dp)                            :: maxlogweight_global

        !> what to thin the posterior by
        real(dp) :: thin_posterior

        !> number of accepted sampled points
        integer :: n_accepted

        !> params of best fit
        real(dp), allocatable, dimension(:) :: params_best_fit

        !> the best loglikelihood so far
        real(dp) :: best_fit

    end type run_time_info

end module
