subroutine SeikInit
    use myalloc
    
    implicit none

    double precision :: TempMask3D(jpk, jpj, jpi), TempMask2D(jpj, jpi)
    integer :: indexi, indexj, indexk
    
    CutLeft=(.not.(nldi==1))
    CutRight=(.not.(nlei==jpi))
    CutTop=(.not.(nlej==jpj))
    CutBottom=(.not.(nldj==1))
    BaseMember=reshape(ModelErrorDiag1,(/jpk,jpj,jpi,jptra/))
    call CutCellsTracer(BaseMember)
    ModelErrorDiag1=reshape(BaseMember,(/SpaceDim/))
    BaseMember=Huge(BaseMember(1,1,1,1)) 

    call readnc_slice_double_2d("BC/TIN_yyyy0115-00:00:00.nc", "riv_N3n", TempMask2D)
    call readnc_slice_double("BC/GIB_yyyy0215-12:00:00.nc", "gib_N6r", TempMask3D)
    SeikMask=0
    do indexi=1,jpi
        do indexj=1, jpj
            do indexk=1, 30 !jpk
                if (abs(TempMask2D(indexj, indexi)+1)>1.0d-6) cycle
                if ((abs(TempMask3D(indexk, indexj, indexi)+1)<1.0d-6).and.(bfmmask(indexk, indexj, indexi)==1)) SeikMask(indexk, indexj, indexi)=1
            end do
        end do
    end do

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
