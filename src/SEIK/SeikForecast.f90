subroutine SeikForecast()
    Use myalloc
    use mpi

    implicit none

    INTEGER :: ierr, indexi

    where (trn<1.0d-12) trn=1.0d-12 !this should already be done with SMALL in another routine, but I don't have time to check now
    do indexi=1, jptra
        where (tmask==1) 
            trnEnsemble(:,:,:,indexi)=trn(:,:,:,indexi)
            BaseMember(:,:,:,indexi)=log(trn(:,:,:,indexi))
            trnEnsembleWeighted(:,:,:,indexi)=BaseMember(:,:,:,indexi)*Weight
        elsewhere
            trnEnsemble(:,:,:,indexi)=0.0d0
            BaseMember(:,:,:,indexi)=0.0d0
            trnEnsembleWeighted(:,:,:,indexi)=0.0d0
        end where
    end do
    !trnEnsemble=log(trn)
    !trnEnsembleWeighted=trnEnsemble*Weight

    call MPI_AllReduce(trnEnsembleWeighted, trn, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    !BaseMember=trnEnsemble-trn
    where (trn>log(1.0d-5)) 
        BaseMember=BaseMember-trn
    elsewhere
        BaseMember=0.0d0
    end where

    !do indexi=1, jptra
        !where (bfmmask==0) BaseMember(:,:,:,indexi)=0.0d0
    !end do

    !if (UseInflation==.true.) then
        !write(*,*) 'Missing code using inflation: check if is it necessary to prepare L and other matrices for future use. I will stop!'
        !call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        !if (GlobalRank==0) CovSeik1=ForgettingFactor*TTTSeik
    !else
    
    if (EnsembleRank==NotWorkingMember) then
        call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        if (MyRank==0) then
            if (UseInflation) then
                CovSeik1=ForgettingFactor*TTTSeik
            else
                call mpi_gatherv(0,0,mpi_real8,LTQ1L,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
                CovSmoother1part=TTTSeik+LTQ1L
                CovSeik1=TTTSeik
                TempMatrixSeik=CovSmoother1part
                call InvMatMul(TempMatrixSeik,CovSeik1,SeikDim,ierr)
                TempMatrixSeik=matmul(LTQ1L,CovSeik1)
                CovSeik1=TempMatrixSeik
            end if
        end if
    else
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        if (.not.UseInflation) then
            TempVecSeik=reshape(BaseMember,(/SpaceDim/))
            TempVecSeik=TempVecSeik*ModelErrorDiag1
            TempSliceSeik=matmul(TempVecSeik,LSeik)
            if (MyRank==0) then
                call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
                call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,0,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
            else
                call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            end if
        end if
    end if
        
    !end if
    
    trn=exp(trn)
    do indexi=1, jptra
        trn(:,:,:,indexi)=trn(:,:,:,indexi)*tmask
    end do
    trb=trn

end subroutine
