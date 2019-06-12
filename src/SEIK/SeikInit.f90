subroutine SeikInit
    use myalloc
    
    implicit none
    
    CutLeft=(.not.(nldi==1))
    CutRight=(.not.(nlei==jpi))
    CutTop=(.not.(nlej==jpj))
    CutBottom=(.not.(nldj==1))
    BaseMember=reshape(ModelErrorDiag1,(/jpk,jpj,jpi,jptra/))
    call CutCellsTracer(BaseMember)
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
    
    double precision, dimension(jpk, jpj, jpi, jptra), intent(inout) :: SurfaceTracer
    
    if (CutLeft) SurfaceTracer(:,1)=0.0d0
    if (CutRight) SurfaceTracer(:,jpi)=0.0d0
    if (CutTop) SurfaceTracer(jpj,:)=0.0d0
    if (CutBottom) SurfaceTracer(1,:)=0.0d0
    
end subroutine
