subroutine ModSeikAnalysis
    use mpi
    use myalloc

    implicit none
    
    integer :: ierr

    ObsBaseMember=ComputedObsSeik*SeikWeight
    call MPI_AllReduce(ObsBaseMember, ComputedObsMean, ObsSpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)

    ObsBaseMember=ComputedObsSeik-ComputedObsMean
    
    !if (CutLeft) ObsBaseMember(:,1)=0.0d0
    !if (CutRight) ObsBaseMember(:,jpi)=0.0d0
    !if (CutTop) ObsBaseMember(jpj,:)=0.0d0
    !if (CutBottom) ObsBaseMember(1,:)=0.0d0
    
    if (EnsembleRank==NotWorkingMember) then
        call mpi_allgatherv(0,0,mpi_real8,HLSeik,MpiCountObs,MpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
        TempObsSeik=reshape(MisfitSeik,(/ObsSpaceDim/))
    else
        call mpi_allgatherv(ObsBaseMember,ObsSpaceDim,mpi_real8,HLSeik,MpiCountObs,MpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
        TempObsSeik=reshape(ObsBaseMember,(/ObsSpaceDim/))
    end if
    
    TempObsSeik=TempObsSeik*ObsErrorDiag1
    TempSliceSeik=matmul(TempObsSeik,HLSeik)
    if (MyRank==0) then
        call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
        if (EnsembleRank==NotWorkingMember) then
            call mpi_gatherv(0,0,mpi_real8,HLTR1HL,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
            CovSeik1=CovSeik1+HLTR1HL
            TempMatrixSeik=CovSeik1
            call dposv( 'U', SeikDim, 1, TempMatrixSeik, SeikDim, TempSliceSeik2, SeikDim, ierr)
            if (ierr.ne.0) error stop 'Analysis Matrix inversion failed!'         
        else
            call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,0,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
        end if
        ChangeCoefSeik=TempSliceSeik2
        call MPI_Bcast( ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)        
    else
        call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
    end if
    
    call MPI_Bcast( ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    TempVecSeik=matmul(Lseik,ChangeCoefSeik)
    BaseMember=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))
    
    if (UseMaxVarSEIK) then
        trnEnsembleWeighted=BaseMember**2
        where (trnEnsembleWeighted>MaxVarSEIK) BaseMember=sign(sqrt(MaxVarSEIK),BaseMember)
    end if
    
    
    BaseMember=exp(BaseMember)
    trn=trn*BaseMember
    !trb=trn
end subroutine
 
