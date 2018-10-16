! iso_c_binding interface for MultiNest_v3.7 Nested::nestRun
! should now be compliant with Fortran 2003 rather than 2008
! Michele Vallisneri, 2014/10/16
!
! allows calling multinest with the "natural" C prototype
!
! void run(bool nest_IS,bool nest_mmodal,bool nest_ceff, \
!          index_t nest_nlive,value_t nest_tol,value_t nest_ef,index_t nest_ndims,index_t nest_totPar,index_t nest_nCdims,index_t maxClst, \
!          index_t nest_updInt,value_t nest_Ztol,char nest_root[],index_t seed,index_t nest_pWrap[], \
!          bool nest_fb,bool nest_resume,bool nest_outfile,bool initMPI,value_t nest_logZero,index_t nest_maxIter, \
!          value_t (*loglike)(value_t *,index_t,index_t,void *context), \
!          void (*dumper)(index_t,index_t,index_t,value_t *,value_t *,value_t *,value_t,value_t,value_t,void *context),void *context);
!
! as well as
!
! value_t loglike(value_t *Cube,index_t n_dim,index_t nPar,void *context);
! void dumper(index_t nSamples,index_t nlive,index_t nPar, \
!             value_t *physLive,value_t *posterior,value_t *paramConstr, \
!             value_t maxLogLike,value_t logZ,value_t INSlogZ,value_t logZerr,void *context);
!
! note that we are assuming that (void *) is the same size as index_t, but that's what multinest uses

module cnested
	use iso_c_binding, only: c_funptr
   	type(c_funptr) :: theloglike, thedumper

	contains

 	subroutine loglike_f(Cube,n_dim,nPar,lnew,context_pass)
	use iso_c_binding, only: c_double, c_f_procpointer

	implicit none

	integer          :: n_dim,nPar,context_pass
	value_t precision :: Cube(nPar)
	value_t precision :: lnew

	interface
		real(c_double) function loglike_proto(Cube,n_dim,nPar,context)
		use iso_c_binding, only: c_int, c_double, c_ptr

		implicit none

		integer(c_int), intent(in), value :: n_dim,nPar
		real(c_double), intent(inout)     :: Cube(nPar)
		integer(c_int), intent(in) :: context
		! better, but "transfer" is problematic:
		! type(c_ptr),  intent(in) :: context

		end function loglike_proto
	end interface

	procedure(loglike_proto), pointer :: loglike_c
	call c_f_procpointer(theloglike,loglike_c)

	! type(c_ptr) :: context_c
	! context_c = transfer(context_pass,context_c)

	lnew = loglike_c(Cube,n_dim,nPar,context_pass)

	end subroutine loglike_f


	subroutine dumper_f(nSamples,nlive,nPar,physLive,posterior,&
        paramConstr,maxLogLike,logZ,INSlogZ,logZerr,n_accepted,context_pass)
	use iso_c_binding, only: c_double, c_f_procpointer

	implicit none

	integer          :: nSamples, nlive, nPar, context_pass, n_accepted
	value_t precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
	value_t precision :: maxLogLike, logZ, INSlogZ, logZerr

	interface
		subroutine dumper_proto(nSamples,nlive,nPar,physLive,posterior,&
            paramConstr,maxLogLike,logZ,INSlogZ,logZerr,n_accepted,context)
		use iso_c_binding, only: c_int, c_double, c_ptr

		implicit none

		integer(c_int), intent(in), value :: nSamples, nlive, nPar, n_accepted
		real(c_double), intent(in) :: physLive(nlive,nPar+1), posterior(nSamples,nPar+2),paramConstr(4*nPar)
		real(c_double), intent(in), value :: maxLogLike, logZ, INSlogZ, logZerr
		integer(c_int), intent(in), value :: context
		! better, but "transfer" is problematic:
		! type(c_ptr),  intent(in) :: context

		end subroutine dumper_proto
	end interface

	procedure(dumper_proto), pointer :: dumper_c
	call c_f_procpointer(thedumper,dumper_c)

	! type(c_ptr) :: context_c
	! context_c = transfer(context_pass,context_c)

	call dumper_c(nSamples,nlive,nPar,physLive,posterior,&
        paramConstr,maxLogLike,logZ,INSlogZ,logZerr,n_accepted,context_pass)

	end subroutine dumper_f


	subroutine run(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar,nest_nCdims,maxClst, &
	nest_updInt,nest_Ztol,nest_root,seed,nest_pWrap,nest_fb,nest_resume,nest_outfile,initMPI,nest_logZero,nest_maxIter, &
	loglike,dumper,context) bind(c)

	use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, c_ptr, C_NULL_CHAR
	use Nested, only: nestRun
	implicit none

	integer(c_int),  intent(in), value :: nest_ndims,nest_nlive,nest_updInt,seed
	integer(c_int),  intent(in), value :: maxClst,nest_totPar,nest_nCdims,nest_maxIter
	integer(c_int),  intent(in) :: nest_pWrap(nest_ndims)
	logical(c_bool), intent(in), value :: nest_IS,nest_mmodal,nest_fb,nest_resume,nest_ceff,nest_outfile,initMPI
	character(kind=c_char,len=1), dimension(1), intent(in) :: nest_root
	real(c_double),  intent(in), value :: nest_tol,nest_ef,nest_Ztol,nest_logZero
   	type(c_funptr),  intent(in), value :: loglike, dumper
	type(c_ptr),     intent(in) :: context

	character(len=100) :: fnest_root
	integer :: i, context_f

	fnest_root = ' '
	do i = 1, 100
		if (nest_root(i) == C_NULL_CHAR) then
			exit
		else
			fnest_root(i:i) = nest_root(i)
		end if
	end do

        theloglike = loglike
        thedumper = dumper

	! context_f = transfer(context,context_f)

	call nestRun(logical(nest_IS),logical(nest_mmodal),logical(nest_ceff), &
	nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar,nest_nCdims,maxClst, &
	nest_updInt,nest_Ztol,fnest_root,seed,nest_pWrap, &
	logical(nest_fb),logical(nest_resume),logical(nest_outfile),logical(initMPI),nest_logZero,nest_maxIter, &
	loglike_f,dumper_f,context_f)

	end subroutine

end module
