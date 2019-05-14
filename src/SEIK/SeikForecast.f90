subroutine SeikForecast()
	Use myalloc
	use mpi

	implicit none

	INTEGER :: ierr

	trnEnsemble=trn
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
                    CovSeik1=matmul(TTTSeik,CovSeik1)
                end if
            else
                call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
                TempVecSeik=reshape(BaseMember,(/SpaceDim/))
                TempVecSeik=TempVecSeik*ModelErrorDiag1
                TempSliceSeik=matmul(TempVecSeik,LSeik)
                call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
                if (MyRank==0) call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,CovSmoother1part,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
            end if
            
  	end if

	trb=trn

end subroutine
