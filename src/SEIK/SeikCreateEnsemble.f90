subroutine SeikCreateEnsemble()
    Use myalloc
    use mpi
    
    implicit none
    integer :: ierr
    
    if (MyRank==0) then
        if (EnsembleRank==NotWorkingMember) call Sampling(CovSeik1, SeikDim, ChangeBaseSeik, ierr)
        call MPI_Scatter(ChangeBaseSeik, SeikDim, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
    end if
    call MPI_Bcast(ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    TempVecSeik=matmul(Lseik,ChangeCoefSeik)
    trn=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))
    trb=trn
end subroutine
