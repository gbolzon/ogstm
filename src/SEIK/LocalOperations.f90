subroutine LocalSendRecive(patch)
    use myalloc
    use mpi
    
    implicit none
    double precision, dimension(SeikDim, jpj, jpi), intent(in) :: patch
    integer :: ierr
    
!               3     4
!               |     ^
!               |     |
!               v     |
!           ________________
!          |                |
!     1<-- |                | 1 <--
!     2--> |                | 2 -->
!          |________________|
!               3     4
!               |     ^
!               |     |
!               v     |    
    
    LocalPatch=0.0d0
    LocalPatch(:,nldj:nlej,nldi:nlei)=patch(:,nldj:nlej,nldi:nlei)

    if (mod(xRank, 2)==0) then
        if (nowe/=-1) then
            HorizontalPatch=patch(:,nldj:nlej,2:2+LocalRange)            
            call MPI_Send(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, nowe, 10, LocalComm, ierr)
            call MPI_Recv(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, nowe, 20, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,nldj:nlej, 1-LocalRange:1)=HorizontalPatch
        end if
        
        if (noea/=-1) then
            HorizontalPatch=patch(:,nldj:nlej,jpi-1-LocalRange:jpi-1)
            call MPI_Send(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, noea, 20, LocalComm, ierr)
            call MPI_Recv(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, noea, 10, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,nldj:nlej, jpi:jpi+LocalRange)=HorizontalPatch
        end if
    else
        if (noea/=-1) then     
            call MPI_Recv(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, noea, 10, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,nldj:nlej, jpi:jpi+LocalRange)=HorizontalPatch
            HorizontalPatch=patch(:,nldj:nlej,jpi-1-LocalRange:jpi-1)
            call MPI_Send(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, noea, 20, LocalComm, ierr)
        end if
    
        if (nowe/=-1) then
            call MPI_Recv(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, nowe, 20, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,nldj:nlej, 1-LocalRange:1)=HorizontalPatch
            HorizontalPatch=patch(:,nldj:nlej,2:2+LocalRange)
            call MPI_Send(HorizontalPatch, (nlej-nldj+1)*(LocalRange+1)*SeikDim, MPI_real8, nowe, 10, LocalComm, ierr)   
        end if
    end if
    
    if (mod(yRank, 2)==0) then
        if (noso/=-1) then
            VerticalPatch=patch(:,2:2+LocalRange,nldi:nlei)
            call MPI_Send(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, noso, 30, LocalComm, ierr)
            call MPI_Recv(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, noso, 40, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,1-LocalRange:1,nldi:nlei)=VerticalPatch
        end if
        
        if (nono/=-1) then
            VerticalPatch=patch(:,jpj-1-LocalRange:jpj-1,nldi:nlei)
            call MPI_Send(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, nono, 40, LocalComm, ierr)
            call MPI_Recv(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, nono, 30, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,jpj:jpj+LocalRange,nldi:nlei)=VerticalPatch
        end if
    else
        if (nono/=-1) then
            call MPI_Recv(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, nono, 30, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,jpj:jpj+LocalRange,nldi:nlei)=VerticalPatch
            VerticalPatch=patch(:,jpj-1-LocalRange:jpj-1,nldi:nlei)
            call MPI_Send(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, nono, 40, LocalComm, ierr)
        end if
    
        if (noso/=-1) then
            call MPI_Recv(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, noso, 40, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,1-LocalRange:1,nldi:nlei)=VerticalPatch
            VerticalPatch=patch(:,2:2+LocalRange,nldi:nlei)
            call MPI_Send(VerticalPatch, (LocalRange+1)*(nlei-nldi+1)*SeikDim, MPI_real8, noso, 30, LocalComm, ierr)
        end if
    end if
    
    if (mod((xRank+yRank)/2, 2)==0) then
        if (nosowe/=-1) then
            DiagonalPatch=patch(:,2:2+LocalRange,2:2+LocalRange)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosowe, 50, LocalComm, ierr)
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosowe, 60, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,1-LocalRange:1, 1-LocalRange:1)=DiagonalPatch
        end if
        
        if (nonoea/=-1) then
            DiagonalPatch=patch(:,jpj-1-LocalRange:jpj-1,jpi-1-LocalRange:jpi-1)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonoea, 60, LocalComm, ierr)
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonoea, 50, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,jpj:jpj+LocalRange, jpi:jpi+LocalRange)=DiagonalPatch
        end if
    else
        if (nonoea/=-1) then
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonoea, 50, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,jpj:jpj+LocalRange, jpi:jpi+LocalRange)=DiagonalPatch
            DiagonalPatch=patch(:,jpj-1-LocalRange:jpj-1,jpi-1-LocalRange:jpi-1)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonoea, 60, LocalComm, ierr)
        end if
    
        if (nosowe/=-1) then
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosowe, 60, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,1-LocalRange:1, 1-LocalRange:1)=DiagonalPatch
            DiagonalPatch=patch(:,2:2+LocalRange,2:2+LocalRange)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosowe, 50, LocalComm, ierr)
        end if
    end if
    
    if (mod(floor((xRank-yRank)/2.0), 2)==0) then
        if (nosoea/=-1) then
            DiagonalPatch=patch(:,2:2+LocalRange,jpi-1-LocalRange:jpi-1)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosoea, 70, LocalComm, ierr)
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosoea, 80, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,1-LocalRange:1, jpi:jpi+LocalRange)=DiagonalPatch
        end if
        
        if (nonowe/=-1) then
            DiagonalPatch=patch(:,jpj-1-LocalRange:jpj-1,2:2+LocalRange)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonowe, 80, LocalComm, ierr)
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonowe, 70, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,jpj:jpj+LocalRange, 1-LocalRange:1)=DiagonalPatch
        end if
    else
        if (nonowe/=-1) then
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonowe, 70, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,jpj:jpj+LocalRange, 1-LocalRange:1)=DiagonalPatch
            DiagonalPatch=patch(:,jpj-1-LocalRange:jpj-1,2:2+LocalRange)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nonowe, 80, LocalComm, ierr)
        end if
    
        if (nosoea/=-1) then
            call MPI_Recv(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosoea, 80, LocalComm, MPI_STATUS_IGNORE, ierr)
            LocalPatch(:,1-LocalRange:1, jpi:jpi+LocalRange)=DiagonalPatch
            DiagonalPatch=patch(:,2:2+LocalRange,jpi-1-LocalRange:jpi-1)
            call MPI_Send(DiagonalPatch, (LocalRange+1)*(LocalRange+1)*SeikDim, MPI_real8, nosoea, 70, LocalComm, ierr)
        end if
    end if
   
end subroutine

subroutine SummingLocalPatch(patch)
    use myalloc
    
    implicit none
    double precision, dimension(SeikDim, jpj, jpi), intent(out) :: patch
    integer :: indexi, indexj, indexk, temp
    
    patch=0.0d0
    do indexi=-LocalRange,LocalRange
        temp=floor(sqrt((0.5d0+LocalRange)**2-indexi**2))
        do indexj=-temp+1, temp+1
            patch(:,1,1)=patch(:,1,1)+LocalPatch(:,indexj, 1+indexi)
        end do
    end do
    do indexj=2, jpj
        patch(:,indexj,1)=patch(:,indexj-1,1)
        do indexi=-LocalRange,LocalRange
            temp=floor(sqrt((0.5d0+LocalRange)**2-indexi**2))
            patch(:,indexj,1)=patch(:,indexj,1)-LocalPatch(:,indexj-1-temp, 1+indexi)+LocalPatch(:,indexj+temp, 1+indexi)
        end do
    end do
    do indexi=2, jpi
        do indexj=1,jpj
            patch(:,indexj,indexi)=patch(:,indexj,indexi-1)
            do indexk=-LocalRange,LocalRange
                temp=floor(sqrt((0.5d0+LocalRange)**2-indexk**2))
                patch(:,indexj,indexi)=patch(:,indexj,indexi)-LocalPatch(:,indexj+indexk, indexi-1-temp)+LocalPatch(:,indexj+indexk, indexi+temp)
            end do
        end do
    end do
    
end subroutine

subroutine SummingLocalPatchWeighted(patch)
    use myalloc
    
    implicit none
    double precision, dimension(SeikDim, jpj, jpi), intent(out) :: patch
    integer :: indexi, indexj, indexk, indexl, temp
    double precision :: localweight
    
    patch=0.0d0
    do indexi=1, jpi
        do indexj=1,jpj
            do indexk=-LocalRange,LocalRange
                temp=floor(sqrt((0.5d0+LocalRange)**2-indexk**2))
                do indexl=-temp,temp
                    call LocalWeightFunction(dble(indexk**2+ indexl**2), localweight)
                    patch(:,indexj,indexi)=patch(:,indexj,indexi)+LocalPatch(:,indexj+indexl,indexi+indexk)*localweight
                end do
            end do
        end do
    end do
    
end subroutine

subroutine LocalProduct(ChangeBase_sji, output)
    use myalloc
    
    implicit none
    double precision, dimension(SeikDim, jpj, jpi), intent(in) :: ChangeBase_sji
    double precision, dimension(jpk, jpj, jpi, jptra), intent(out) :: output
    integer :: indexi, indexj, indexk
    
    do indexi=1, jpi
        do indexj=1, jpj
            if (SeikMask(1,indexj, indexi)==1) then
                output(:, indexj, indexi, :)=ChangeBase_sji(1,indexj, indexi)*LSeik_reshape(:,indexj, indexi, :, 1)
                do indexk=2, SeikDim
                    output(:, indexj, indexi, :)=output(:, indexj, indexi, :)+ChangeBase_sji(indexk,indexj, indexi)*LSeik_reshape(:,indexj, indexi, :, indexk)
                end do
            else
                output(:, indexj, indexi, :)=0.0d0
            end if
        end do
    end do
    
end subroutine
