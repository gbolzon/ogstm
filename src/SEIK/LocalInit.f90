subroutine LocalInit(domdec, xpos, ypos)
    Use myalloc

    implicit none
    integer, dimension(CommSize,13), intent(in) :: domdec
    integer, intent(in) :: xpos, ypos
    
    xRank=xpos
    yRank=ypos
    
    if (noea==-1) then
        if (nono==-1) then 
            nonoea=-1
        else
            nonoea=domdec(nono+1, 11)
        end if
        
        if (noso==-1) then 
            nosoea=-1
        else
            nosoea=domdec(noso+1, 11)
        end if
    else
        nonoea=domdec(noea+1, 12)
        nosoea=domdec(noea+1, 13)
    end if
    
    if (nowe==-1) then
        if (nono==-1) then 
            nonowe=-1
        else
            nonowe=domdec(nono+1, 10)
        end if
        
        if (noso==-1) then 
            nosowe=-1
        else
            nosowe=domdec(noso+1, 10)
        end if
    else
        nonowe=domdec(nowe+1, 12)
        nosowe=domdec(nowe+1, 13)
    end if
    
    LocalMaxRange=minval(domdec(:,4:5))-3
    
end subroutine
