subroutine Sampling(CovMatrix1, nside, ChangeBase, ierr)
    use mpi
    use myalloc
    
    implicit none
    
    integer, intent(in) :: nside
    integer, intent(out) :: ierr
    double precision, dimension(nside,nside), intent(inout) :: CovMatrix1
    double precision, dimension(nside,0:nside), intent(out) :: ChangeBase
    integer :: indexi
    
    call dpotrf( 'U', nside, CovMatrix1, nside, ierr )
    if (ierr.ne.0) error stop 'Choleky failed'
    
    OrtMatrixSampling(:,1)=AllWeightsSqrt
    call OrtMatrix(OrtMatrixSampling, SeikDim+1, 1, SeikDim)
    do indexi=2, SeikDim+1
        OrtMatrixSampling(:,indexi)=OrtMatrixSampling(:,indexi)*AllWeightsSqrt1
    end do
    ChangeBase=transpose(OrtMatrixSampling(:,2:SeikDim+1))
    call dtrtrs( 'U', 'N', 'N', SeikDim, SeikDim+1, CovMatrix1, SeikDim, ChangeBase, SeikDim, ierr)
end subroutine
