subroutine SeikCreateEnsemble()
    Use myalloc
    use mpi
    
    implicit none
    integer :: ierr
    
    if (UseModSeik) then
    
    end if
    
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

    trnEnsembleWeighted=BaseMember**2
    trnEnsembleWeighted=trnEnsembleWeighted*Weight
    call MPI_AllReduce(trnEnsembleWeighted, trnVariance, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    where (trnVariance>MaxVarSEIK)
        trnEnsemble=sqrt(MaxVarSEIK/trnEnsemble)
    elsewhere
        trnEnsemble=1.0d0
    end where
    BaseMember=BaseMember*trnEnsemble

    BaseMember=exp(BaseMember)
    trn=trn*BaseMember
    trb=trn
end subroutine
