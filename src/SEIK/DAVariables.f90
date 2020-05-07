module DAVariables
    implicit none
    
    contains
    
    function IsHighVariance(name)
        implicit none
        
        character(LEN=*), intent(in) :: name
        logical :: IsHighVariance
        
        if (name(1:1).eq."N") then
            IsHighVariance=.false.
        else
            IsHighVariance=.true.
        end if
    end function
    
end module
