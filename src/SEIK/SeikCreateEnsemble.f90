subroutine SeikCreateEnsemble()
    Use myalloc
    use mpi
    
    implicit none
    integer :: ierr
    
    if (MyRank==0) then
        if (EnsembleRank==NotWorkingMember) then
            call Sampling(CovSeik1, SeikDim, ChangeBaseSeik, ierr)
            call MPI_Scatter(ChangeBaseSeik, SeikDim, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        else
            call MPI_Scatter(0, 0, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        end if
    end if
    call MPI_Bcast(ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    TempVecSeik=matmul(Lseik,ChangeCoefSeik)
    BaseMember=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))
    BaseMember=exp(BaseMember)
    trn=trn*BaseMember
    trb=trn
end subroutine
