subroutine AnalysisWrapper(DateString)
    use myalloc
    use mpi
    
    implicit none
    character(LEN=17) :: DateString
    integer :: indexi, ierr
    
    if ((UseLocalForecast.or.FirstTimeSampling).and.(SeikDim>0)) then
        trb=trn !forse non serve
        call SeikCreateEnsemble
        trnEnsemble=trn
        trn=trb
        
        if (EnsembleRank==NotWorkingMember) then
            call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
            CovSeik1=TTTSeik
        else
            where (SeikTrcMask==1)
                BaseMember=log(trnEnsemble)-log(trn)
            elsewhere
                BaseMember=0.0d0
            end where
            call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        end if
        
    end if
    
    CALL MainAssimilationSeik(DateString) 

end subroutine
