subroutine SeikAnalysis
    use mpi
    use myalloc

    implicit none
    
    integer :: ierr

    ObsBaseMember=ComputedObsSeik*Weight
    call MPI_AllReduce(ObsBaseMember, ComputedObsMean, ObsSpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    
    if (EnsembleRank==NotWorkingMember) then
        ObsBaseMember=MisfitSeik
        call mpi_allgatherv(0,0,mpi_real8,HLSeik,MpiCountObs,MpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
    else
        ObsBaseMember=ComputedObsSeik-ComputedObsMean
    
        !if (CutLeft) ObsBaseMember(:,1)=0.0d0
        !if (CutRight) ObsBaseMember(:,jpi)=0.0d0
        !if (CutTop) ObsBaseMember(jpj,:)=0.0d0
        !if (CutBottom) ObsBaseMember(1,:)=0.0d0
        
        call mpi_allgatherv(ObsBaseMember,ObsSpaceDim,mpi_real8,HLSeik,MpiCountObs,MpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
    end if
    
    TempObsSeik=reshape(ObsBaseMember,(/ObsSpaceDim/))
    TempObsSeik=TempObsSeik*ObsErrorDiag1
    TempSliceSeik=matmul(TempObsSeik,HLSeik)
    
    if (UseDiffCov) then
        UDiffObsBaseMember=(ObsBaseMember(:,2:jpi)-ObsBaseMember(:,1:jpi-1))/e1u(:,1:jpi-1)*SeikUMask(1,:,:)
        VDiffObsBaseMember=(ObsBaseMember(2:jpi,:)-ObsBaseMember(1:jpi-1,:))/e2v(1:jpj-1,:)*SeikVMask(1,:,:)
        
        if (EnsembleRank==NotWorkingMember) then
            call mpi_allgatherv(0,0,mpi_real8,UHLSeik,UDiffMpiCountObs,UDiffMpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
            call mpi_allgatherv(0,0,mpi_real8,VHLSeik,VDiffMpiCountObs,VDiffMpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
        else            
            call mpi_allgatherv(UDiffObsBaseMember,UDiffObsSpaceDim,mpi_real8,UHLSeik,UDiffMpiCountObs,UDiffMpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
            call mpi_allgatherv(VDiffObsBaseMember,VDiffObsSpaceDim,mpi_real8,VHLSeik,VDiffMpiCountObs,VDiffMpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
        end if
        
        TempUDiffObsSeik=reshape(UDiffObsBaseMember,(/UDiffObsSpaceDim/))
        TempUDiffObsSeik=TempUDiffObsSeik*UDiffObsErrorDiag1
        TempSliceSeik=TempSliceSeik+matmul(TempUDiffObsSeik,UHLSeik)
        
        TempVDiffObsSeik=reshape(VDiffObsBaseMember,(/VDiffObsSpaceDim/))
        TempVDiffObsSeik=TempVDiffObsSeik*VDiffObsErrorDiag1
        TempSliceSeik=TempSliceSeik+matmul(TempVDiffObsSeik,VHLSeik)
    end if
    
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
        ChangeCoefSeik=TempSliceSeik2 !Si puo' eliminare TempSlice2 in facore di ChangeCoefSeik?
        call MPI_Bcast( ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)        
    else
        call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
    end if
    
    call MPI_Bcast( ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    TempVecSeik=matmul(Lseik,ChangeCoefSeik)
    trnEnsemble=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))
    
    if (UseMaxVarSEIK) then
        trnEnsembleWeighted=trnEnsemble**2
        where (trnEnsembleWeighted>MaxVarSEIK) trnEnsemble=sign(sqrt(MaxVarSEIK),trnEnsemble)
    end if
    
    trnEnsemble=exp(trnEnsemble)
    trn=trn*trnEnsemble
    !trb=trn
end subroutine
 