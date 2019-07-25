subroutine SeikInit
    use myalloc
use mpi
    
    implicit none
integer ierr

    double precision :: TempMask3D(jpk, jpj, jpi), TempMask2D(jpj, jpi)
    integer :: indexi, indexj, indexk, indexip1, indexim1, indexjp1, indexjm1
    
    CutLeft=(.not.(nldi==1))
    CutRight=(.not.(nlei==jpi))
    CutTop=(.not.(nlej==jpj))
    CutBottom=(.not.(nldj==1))
    BaseMember=reshape(ModelErrorDiag1,(/jpk,jpj,jpi,jptra/))
    call CutCellsTracer(BaseMember)
    
    call readnc_slice_double_2d("BC/TIN_yyyy0115-00:00:00.nc", "riv_N3n", TempMask2D)
    call readnc_slice_double("BC/GIB_yyyy0215-12:00:00.nc", "gib_N6r", TempMask3D)
    SeikMask=0
    do indexi=1,jpi
        if (indexi==1) then
            indexim1=1
        else
            indexim1=indexi-1
        end if
        if (indexi==jpi) then
            indexip1=jpi
        else
            indexip1=indexi+1
        end if
        do indexj=1, jpj
            if (indexj==1) then
                indexjm1=1
            else
                indexjm1=indexj-1
            end if
            if (indexj==jpj) then
                indexjp1=jpj
            else
                indexjp1=indexj+1
            end if
            do indexk=1, 30 !jpk-1
                if ((abs(TempMask2D(indexj, indexi)+1)>1.0d-6).and.(indexk==1)) cycle
                if (tmask(indexk+1, indexj, indexi)+tmask(indexk, indexjp1, indexi)+tmask(indexk, indexjm1, indexi) & 
                    +tmask(indexk, indexj, indexip1)+tmask(indexk, indexj, indexim1)<5) cycle
                if ((abs(TempMask3D(indexk, indexj, indexi)+1)<1.0d-6).and.(bfmmask(indexk, indexj, indexi)==1)) SeikMask(indexk, indexj, indexi)=1
            end do
        end do
    end do

call  mpi_barrier(mpi_comm_World,ierr)
write (*,*) "e=", EnsembleRank, "m=", MyRank, "t=", sum(tmask), "b=", sum(bfmmask), "s=", sum(SeikMask)
!call  mpi_barrier(mpi_comm_World,ierr)
!if ((EnsembleRank==0).and.(MyRank==27)) then
!write (*,*) "t="
!write (*,*) tmask
!write (*,*) "b="
!write (*,*) bfmmask
!write (*,*) "s="
!write (*,*) SeikMask
!end if
!call  mpi_barrier(mpi_comm_World,ierr)
!stop

    do indexi=1, jptra
        BaseMember(:,:,:,indexi)=BaseMember(:,:,:,indexi)*SeikMask
    end do

!call trcwriSeik("12345678901234567",-1,"RESTARTS/",BaseMember)

    ModelErrorDiag1=reshape(BaseMember,(/SpaceDim/))
    BaseMember=Huge(BaseMember(1,1,1,1))


end subroutine

subroutine CutCellsTracer(FullTracer)
    use myalloc
    
    implicit none
    
    double precision, dimension(jpk, jpj, jpi, jptra), intent(inout) :: FullTracer
    
    if (CutLeft) FullTracer(:,:,1,:)=0.0d0
    if (CutRight) FullTracer(:,:,jpi,:)=0.0d0
    if (CutTop) FullTracer(:,jpj,:,:)=0.0d0
    if (CutBottom) FullTracer(:,1,:,:)=0.0d0
    
end subroutine

subroutine CutCellsSurface(SurfaceTracer)
    use myalloc
    
    implicit none
    
    double precision, dimension(jpj, jpi), intent(inout) :: SurfaceTracer
    
    if (CutLeft) SurfaceTracer(:,1)=0.0d0
    if (CutRight) SurfaceTracer(:,jpi)=0.0d0
    if (CutTop) SurfaceTracer(jpj,:)=0.0d0
    if (CutBottom) SurfaceTracer(1,:)=0.0d0
    
end subroutine
