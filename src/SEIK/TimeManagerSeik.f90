module TimeManagerSeik
    use Time_manager
    Use mpi
    implicit none
    
    type(dump_container) :: ForecastTimes
    
contains

    subroutine LoadTimeManagerSeik
        implicit none
        
        ForecastTimes%FileName = 'forecastTimes'
        ForecastTimes%Name='...'
        call Load_Dump_container(ForecastTimes)
        
        call ForecastSubsetCheck(RESTARTS)
        call ForecastSubsetCheck(da_Times)
        
    end subroutine
    
    subroutine UnLoadTimeManagerSeik
        implicit none
        
        call unload_Dump_container(ForecastTimes)
    end subroutine
        
    LOGICAL FUNCTION IsaForecast(datestring)
        IMPLICIT NONE
        
        CHARACTER(LEN=17), INTENT(IN) :: datestring
        ! LOCAL
        INTEGER I

        IsaForecast = .false.
        
        if (datestring.eq.DATESTART) then
            IsaForecast = .true.
            return
        end if
        
        DO I=1, RESTARTS%N
            if (datestring.eq.RESTARTS%TimeStrings(I)) then
                IsaForecast = .true.
                return
            endif
        ENDDO

    END FUNCTION IsaForecast
    
    subroutine ForecastSubsetCheck(struct)
        IMPLICIT NONE
        
        type(dump_container), INTENT(IN) :: struct
        ! LOCAL
        INTEGER indexi, ierr
        
        DO indexi=1, struct%N
            if ((struct%TimeStrings(indexi).lt.DATESTART).or.(struct%TimeStrings(indexi).gt.DATE__END)) return
            if (.not.isaforecast(struct%TimeStrings(indexi))) then
                write(*,*) "Error! ", struct%TimeStrings(indexi), " is in ", struct%FileName, " but it is not a forecast time!"
                call mpi_abort(mpi_comm_world,1,ierr)
            endif
        ENDDO


    END subroutine 
        
end module
