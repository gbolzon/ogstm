subroutine SeikForecast()
    Use myalloc
    use mpi

    implicit none

    INTEGER :: ierr, indexi

    where (trn<1.0d-11) trn=1.0d-11
    do indexi=1, jptra
        where (tmask==1) 
           trnEnsemble(:,:,:,indexi)=log(trn(:,:,:,indexi))
        elsewhere
            trnEnsemble(:,:,:,indexi)=0.0d0
        end where
    end do
    !trnEnsemble=log(trn)
    trnEnsembleWeighted=trnEnsemble*Weight

    call MPI_AllReduce(trnEnsembleWeighted, trn, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    BaseMember=trnEnsemble-trn

    if (UseInflation==.true.) then
        write(*,*) 'Missing code using inflation: check if is it necessary to prepare L and other matrices for future use. I will stop!'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        if (GlobalRank==0) CovSeik1=ForgettingFactor*TTTSeik
    else
        if (EnsembleRank==NotWorkingMember) then
            call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
            if (MyRank==0) then
                call mpi_gatherv(0,0,mpi_real8,LTQ1L,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
                CovSmoother1part=TTTSeik+LTQ1L
                CovSeik1=LTQ1L
                TempMatrixSeik=CovSmoother1part
                call InvMatMul(TempMatrixSeik,CovSeik1,SeikDim,ierr)
                TempMatrixSeik=matmul(TTTSeik,CovSeik1)
                CovSeik1=TempMatrixSeik
            end if
        else
            call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
            TempVecSeik=reshape(BaseMember,(/SpaceDim/))
            TempVecSeik=TempVecSeik*ModelErrorDiag1
            TempSliceSeik=matmul(TempVecSeik,LSeik)
            if (MyRank==0) then
                call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            else
                call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            end if
            if (MyRank==0) call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,0,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
        end if
        
    end if
    
    trn=exp(trn)
    do indexi=1, jptra
        trn(:,:,:,indexi)=trn(:,:,:,indexi)*tmask
    end do
    trb=trn

end subroutine
