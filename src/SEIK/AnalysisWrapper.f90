subroutine AnalysisWrapper(DateString)
    use myalloc
    
    implicit none
    character(LEN=17) :: DateString
    integer :: indexi
    
    if ((UseLocalForecast.or.FirstTimeSampling).and.(SeikDim>0)) then
        trb=trn !forse non serve
        call SeikCreateEnsemble
        trnEnsemble=trn
        trn=trb
        
        if (EnsembleRank==NotWorkingMember) then
            CovSeik1=0.0d0
            do indexi=1, SeikDim
                CovSeik1(indexi,indexi)=1.0d0
            end do
        end if
        
    end if
    
    CALL MainAssimilationSeik(DateString) 

end subroutine
