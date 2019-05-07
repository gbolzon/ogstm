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
            if (GlobalRank==0) CovSeik1=ForgettingFactor*TTTSeik
  	else
            write(*,*) 'Missing code to use model error instead of inflation. I will stop!'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            
            if (EnsembleRank==NotWorkingMember) then
                call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
            else
                call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
            end if
            if (EnsembleRank /= NotWorkingMember) then
                TempVecSeik=reshape(BaseMember,(/SpaceDim/))
                TempVecSeik=TempVecSeik*ModelErrorDiag1
                TempSliceSeik=matmul(TempVecSeik,LSeik)
                call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
                if (MyRank==0) call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,CovSmoother1part,MpiCountSeik,MpiDisplacementSeik,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
            else
                if (MyRank==0) then
                    call mpi_gatherv(0,0,mpi_real8,CovSmoother1part,MpiCountSeik,MpiDisplacementSeik,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
                    CovSmoother1part=TTT+CovSmoother1part
                
                end if
            end if
            
            !LTQ1Lpart
            
  	end if
  	
  	
		!devi fare la covarianza
	  !call MPI_Reduce((trnEnsemble-trn)**2(*WeightSeik), trn, size(trnEnsemble), mpi_real8, MPI_SUM, 0, EnsembleComm,ierr)

	trb=trn
          !ricorda poi che trb(:,:,:,jn) = trn(:,:,:,jn) uno dei due moltiplicato per tmask

end subroutine
