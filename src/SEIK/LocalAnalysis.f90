subroutine LocalAnalysis
    use mpi
    use myalloc

    implicit none
    
    integer :: ierr, indexi, indexj
integer indexk

    ObsBaseMember=ComputedObsSeik*SeikWeight
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
    
    do indexi=1, SeikDim
        TempObsSeik=reshape(ObsBaseMember,(/ObsSpaceDim/))
        TempObsSeik=TempObsSeik*ObsErrorDiag1
        TempObsSeik=TempObsSeik*HLSeik(:,indexi)
        BaseMember_sji(indexi,:,:)=reshape(TempObsSeik,(/jpj,jpi/))
    end do
        
    call LocalSendRecive(BaseMember_sji)
    
    if (UseLocalObsDumping) then
        call SummingLocalPatchWeighted(BaseMember_sji)
    else
        call SummingLocalPatch(BaseMember_sji)
    end if
    
    if (EnsembleRank==NotWorkingMember) then
    
        call mpi_gatherv(0,0,mpi_real8,HLTR1HL_sjis,LocalMpiCountCov,LocalMpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
        

if (myRank==40) then
    write(*,*) "---------------------------------------------------------"
    write(*,*) "HLTR1HL_sjis"
    write(*,*) "rank 40, x=5, y=5"
    do indexi=1, SeikDim
        write(*,*) HLTR1HL_sjis(:,5,5, indexi)
    end do
    write(*,*) "---------------------------------------------------------"
else
    call mpi_barrier(LocalComm, ierr)
end if

        do indexi=1, jpi
            do indexj=1, jpj
                if (SeikMask(1, indexj,indexi)==1) then
                
                    !HLTR1HL=HLTR1HL_sjis(:,indexj,indexi,:)
                    ChangeCoefSeik=BaseMember_sji(:,indexj, indexi)
                    
if ((myRank==40).and.(indexi==5).and.(indexj==5)) then
write(*,*) "---------------------------------------------------------"
write(*,*) "CovSeik1"
write(*,*) "rank 40, x=5, y=5"
do indexk=1, SeikDim
    write(*,*) CovSeik1(:, indexk)
end do
write(*,*) "---------------------------------------------------------"
end if
    
                    !CovSeik1=CovSeik1+HLTR1HL
                    TempMatrixSeik=CovSeik1+HLTR1HL_sjis(:,indexj,indexi,:)
                    
if ((myRank==40).and.(indexi==5).and.(indexj==5)) then
write(*,*) "---------------------------------------------------------"
write(*,*) "CovSeik1post"
write(*,*) "rank 40, x=5, y=5"
do indexk=1, SeikDim
    write(*,*) TempMatrixSeik(:, indexk) !CovSeik1
end do
write(*,*) "---------------------------------------------------------"
write(*,*) "ChangeCoefSeik"
write(*,*) "rank 40, x=5, y=5"
write(*,*) ChangeCoefSeik
write(*,*) "---------------------------------------------------------"
end if

if (.false.) then
                    TempMatrixSeik=CovSeik1
                    call dposv( 'U', SeikDim, 1, TempMatrixSeik, SeikDim, ChangeCoefSeik, SeikDim, ierr)
                    if (ierr.ne.0) then
                        write(*,*) myRank, indexi, indexj, ': Analysis Matrix inversion failed! Error ', ierr  

do indexk=1, SeikDim
    write(*,*) myRank,":", CovSeik1(:, indexk)
end do
                        error stop 'Analysis Matrix inversion failed!'         
                    end if
                    
                    BaseMember_sji(:,indexj, indexi)=ChangeCoefSeik
                    
                    call SymChangeBase(CovSeik1,indexi, indexj)
                        
                    HLTR1HL_sjis(:,indexj,indexi,:)=CovSeik1
end if

                    call SymChangeBase(TempMatrixSeik,indexi, indexj)
                    
                    BaseMember_sji(:,indexj, indexi)=matmul(TempMatrixSeik,ChangeCoefSeik)
                    HLTR1HL_sjis(:,indexj,indexi,:)=TempMatrixSeik
                    
                else
                    BaseMember_sji(:,indexj, indexi)=0.0d0
                    HLTR1HL_sjis(:,indexj,indexi,:)=0.0d0
                end if
            end do
        end do
        
if (myRank==40) call mpi_barrier(LocalComm, ierr)

        call MPI_Scatterv(HLTR1HL_sjis, LocalMpiCountCov, LocalMpiDisplacementCov, mpi_real8, 0, 0, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        
        CovSeik1=0.0d0
        do indexi=1,SeikDim
            CovSeik1(indexi,indexi)=1.0d0
        end do
        
        call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)

    else
    
        call mpi_gatherv(BaseMember_sji,SeikDim*jpj*jpi,mpi_real8,0,LocalMpiCountCov,LocalMpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
        
        call MPI_Scatterv(0, LocalMpiCountCov, LocalMpiDisplacementCov, mpi_real8, BaseMember_sji, SeikDim*jpj*jpi, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
            
        call LocalProduct(BaseMember_sji, BaseMember)
        
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        
    end if
    
    call MPI_Bcast(BaseMember_sji, SeikDim*jpj*jpi, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
    
    call LocalProduct(BaseMember_sji, trnEnsemble)
    
    if (UseMaxVarSEIK) then
        trnEnsembleWeighted=trnEnsemble**2
        where (trnEnsembleWeighted>MaxVarVec) trnEnsemble=sign(sqrt(MaxVarVec),trnEnsemble)
    end if
    
    trnEnsemble=exp(trnEnsemble)
    trn=trn*trnEnsemble
    !trb=trn
        
if (.false.) then
    if (EnsembleRank==NotWorkingMember) then
        
        call MPI_Scatterv(HLTR1HL_sjis, LocalMpiCountCov, LocalMpiDisplacementCov, mpi_real8, 0, 0, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        
        CovSeik1=0.0d0
        do indexi=1,SeikDim
            CovSeik1(indexi,indexi)=1.0d0
        end do
        
        call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        
    else
        
        call MPI_Scatterv(0, LocalMpiCountCov, LocalMpiDisplacementCov, mpi_real8, BaseMember_sji, SeikDim*jpj*jpi, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
            
        call LocalProduct(BaseMember_sji, BaseMember)
        
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        
    end if
end if
end subroutine
 
