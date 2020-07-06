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
