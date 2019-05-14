subroutine SeikCreateEnsemble()
    Use myalloc
    use mpi
    
    implicit none
    integer :: ierr
    
    if (MyRank==0) then
        if (EnsembleRank==NotWorkingMember) call Sampling(CovSeik1, SeikDim, ChangeBaseSeik, ierr)
        MPI_Scatter(ChangeBaseSeik, SeikDim, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
    end if
    MPI_Bcast(ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    trn=matmul(Lseik,ChangeCoefSeik)
    trb=trn
end subroutine
