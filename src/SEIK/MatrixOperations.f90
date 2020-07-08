subroutine ForecastMatrixOp
    use myalloc
    use mpi
    implicit none
    integer :: ierr, indexi
    
    CovSmoother1part=TTTSeik+LTQ1L
    CovSeik1=TTTSeik
    TempMatrixSeik=CovSmoother1part
    call InvMatMul(TempMatrixSeik,CovSeik1,SeikDim,ierr)
    
    if (ierr/=0) then
        write(*,*) "myRank=", myrank, "EnsembleRank=", EnsembleRank, "inversion failed of TTT+LTQ1L. TTT:"
        do indexi=1, SeikDim
            write(*,*) myrank, ",", EnsembleRank, ":", TTTSeik(:,indexi)
        end do
        write(*,*) myrank, ",", EnsembleRank, ", LTQ1L:"
        do indexi=1, SeikDim
            write(*,*) myrank, ",", EnsembleRank, ":", LTQ1L(:,indexi)
        end do
        write(*,*) myrank, ",", EnsembleRank, ", CovSmoother1part:"
        do indexi=1, SeikDim
            write(*,*) myrank, ",", EnsembleRank, ":", CovSmoother1part(:,indexi)
        end do
        call mpi_abort(Mpi_comm_world,1, ierr)
    end if
    
    TempMatrixSeik=matmul(LTQ1L,CovSeik1)
    CovSeik1=TempMatrixSeik
    
end subroutine

subroutine SymChangeBase(matrix, positionx, positiony)
    Use myalloc
    use mpi
    
    implicit none
    double precision, dimension(SeikDim, SeikDim), intent(inout) :: matrix
integer, intent(in) :: positionx, positiony
    integer :: ierr, indexi, neigenvalues
    double precision dlamch
    
!TempMatrixSeik=matrix !usare tempmatrixseik qui e' rischioso, se non te ne accorgi e la usi fuori dalla procedura diventi matto
    
    call dsyevr("V", "A", "U", SeikDim, matrix, SeikDim, 0.0d0, 0.0d0,0.0d0, 0.0d0, &
        dlamch('S'), neigenvalues, eigenvalues, eigenvectors, SeikDim, &
        isuppz, work, lwork, iwork, liwork, ierr)

    if (ierr/=0) then
        write(*,*) "something wrong with svd. I will stop"
        call mpi_abort(mpi_comm_world,1,ierr)
    end if

    if (SeikDim/=neigenvalues) then
        write(*,*) "something strange in the number of eigenvalues!"
    end if
    
    if (eigenvalues(1)<10**-10) then
        write(*,*) myrank, ",", EnsembleRank,",",positionx,",", positiony ,", SymChangeBase small eigenvalue:"
        write(*,*) myrank, ",", EnsembleRank, ":", eigenvalues
!write(*,*) myrank, ",", EnsembleRank, ", CovSeik1:"
!do indexi=1, SeikDim
!    write(*,*) myrank, ",", EnsembleRank, ":", TempMatrixSeik(:,indexi)
!end do
        call mpi_abort(mpi_comm_world, 1, ierr)
    end if
    
    do indexi=1, SeikDim
        matrix(indexi,:)=eigenvectors(:,indexi)/sqrt(eigenvalues(indexi))
    end do
    matrix=MatMul(eigenvectors,matrix)

end subroutine
