subroutine SeikForecast()
	Use myalloc
	use mpi

	implicit none

	INTEGER :: ierr

	trnEnsemble=trn
  	call MPI_AllReduce(trnEnsemble*WeightSeik, trn, size(trnEnsemble), mpi_real8, MPI_SUM, EnsembleComm,ierr)
		!devi fare la covarianza
	  !call MPI_Reduce((trnEnsemble-trn)**2(*WeightSeik), trn, size(trnEnsemble), mpi_real8, MPI_SUM, 0, EnsembleComm,ierr)
	trb=trn
          !ricorda poi che trb(:,:,:,jn) = trn(:,:,:,jn) uno dei due moltiplicato per tmask

end subroutine
