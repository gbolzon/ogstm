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
    
    function IsSEIKVariable(name)
        implicit none
        
        character(LEN=*), intent(in) :: name
        logical :: IsSEIKVariable
        
        if ((name(1:1).eq."P").or.(name(1:1).eq."N").or.(name(1:1).eq."O")) then
            IsSEIKVariable=.true.
        else
            IsSEIKVariable=.false.
        end if
    end function
    
end module
