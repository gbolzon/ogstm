subroutine Sampling(CovMatrix1, nside, ChangeBase, ierr)
    !use mpi
    use myalloc
    
    implicit none
    
    integer, intent(in) :: nside
    integer, intent(out) :: ierr
    double precision, dimension(nside,nside), intent(inout) :: CovMatrix1
    double precision, dimension(nside,0:nside), intent(out) :: ChangeBase
    integer :: indexi
    double precision :: dlamch

    if (UseCholesky) then
        call dpotrf( 'U', nside, CovMatrix1, nside, ierr )
        if (ierr.ne.0) error stop 'Cholesky failed!'
    else
        call dsyevr("V", "A", "U", nside, CovMatrix1, nside, 0.0d0, 0.0d0,0.0d0, 0.0d0, &
            dlamch('S'), neigenvalues, eigenvalues, eigenvectors, nside, &
            isuppz, work, nside*nside, iwork, nside*nside, ierr)

        if (ierr/=0) then
            write(*,*) "something wrong with svd. I will stop"
            call mpi_abort(mpi_comm_world,1,ierr)
        end if

        if (SeikDim/=neigenvalues) then
            write(*,*) "something strange in the number of eigenvalues!"
            write(*,*) "maxNeigenvectors=", maxNeigenvectors, " neigenvalues=", neigenvalues
        end if

        do indexi=1, SeikDim
            CovSeik1(:,indexi)=eigenvectors(:,indexi)/sqrt(eigenvalues(indexi))
        end do
        CovSeik1=Inverse(CovSeik1)
        CovSeik1=MatMul(eigenvectos,CovSeik1)
    end if
    
!!! ATTENZIONE!!! qui ci sono seikdim al posto di nside, verifica!!!!
    
    OrtMatrixSampling(:,1)=AllWeightsSqrt
    call OrtMatrix(OrtMatrixSampling, SeikDim+1,SeikDim+1, 1, SeikDim) 
    do indexi=2, SeikDim+1
        OrtMatrixSampling(:,indexi)=OrtMatrixSampling(:,indexi)*AllWeightsSqrt1
    end do
    ChangeBase=transpose(OrtMatrixSampling(:,2:SeikDim+1))
    

    if (UseCholesky) then
        call dtrtrs( 'U', 'N', 'N', SeikDim, SeikDim+1, CovMatrix1, SeikDim, ChangeBase, SeikDim, ierr)
        if (ierr.ne.0) error stop 'Sampling inversion failed'
    else
        ChangeBase=MatMul(CovMatrix1,ChangeBase)
    end if

end subroutine
