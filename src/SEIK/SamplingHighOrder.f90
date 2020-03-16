subroutine SamplingHighOrder(nside, ChangeBase)
    !use mpi
    use myalloc
    
    implicit none
    
    integer, intent(in) :: nside
    !integer, intent(out) :: ierr
    !double precision, dimension(nside,nside), intent(inout) :: CovMatrix1
    double precision, dimension(nside,0:nside), intent(out) :: ChangeBase
    integer :: indexi
    !double precision :: dlamch

    call OrtMatrix(HighOrderRotation, HighOrderDim,HighOrderDim, 0, HighOrderDim)    
    OrtMatrixSampling(:,1)=HighOrderMatrix(:,1)
    OrtMatrixSampling(:,2:HighOrderDim+1)=MatMul(HighOrderMatrix(:,2:HighOrderDim+1),HighOrderRotation)
    call OrtMatrix(OrtMatrixSampling, nside+1,nside+1, HighOrderDim+1, nside-HighOrderDim) 
    do indexi=2, nside+1
        OrtMatrixSampling(:,indexi)=OrtMatrixSampling(:,indexi)*AllWeightsSqrt1
    end do
    ChangeBase=transpose(OrtMatrixSampling(:,2:nside+1))

end subroutine
