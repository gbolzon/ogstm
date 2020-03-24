subroutine SeikInit
    use myalloc
    use mpi
    USE DA_mem
    
    implicit none
    
    integer ierr
    double precision :: TempMask3D(jpk, jpj, jpi), TempMask2D(jpj, jpi)
    integer :: indexi, indexj, indexk, indexip1, indexim1, indexjp1, indexjm1
    double precision :: partialsum, totalsum

    CutLeft=(.not.(nldi==1))
    CutRight=(.not.(nlei==jpi))
    CutTop=(.not.(nlej==jpj))
    CutBottom=(.not.(nldj==1))

    
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
            do indexk=1, jpk-1
                if ((abs(TempMask2D(indexj, indexi)+1)>1.0d-6).and.(indexk==1)) cycle
                if (tmask(indexk+1, indexj, indexi)+tmask(indexk, indexjp1, indexi)+tmask(indexk, indexjm1, indexi) & 
                    +tmask(indexk, indexj, indexip1)+tmask(indexk, indexj, indexim1)<5) cycle
                if ((abs(TempMask3D(indexk, indexj, indexi)+1)<1.0d-6).and.(bfmmask(indexk, indexj, indexi)==1)) SeikMask(indexk, indexj, indexi)=1
            end do
        end do
    end do

    SeikMask=BfmMask !per il momento preferisco la bfmmask
    SeikTrcMask=0
    do indexi=1,jptra
        !if (isaDAVAR(ctrcnm(indexi))) SeikTrcMask(:,:,:,indexi)=SeikMask usa questa se vuoi solo le variabili di clorofilla
        SeikTrcMask(:,:,:,indexi)=SeikMask
    end do


!call  mpi_barrier(mpi_comm_World,ierr)
!write (*,*) "e=", EnsembleRank, "m=", MyRank, "t=", sum(tmask), "b=", sum(bfmmask), "s=", sum(SeikMask)
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
        do indexj=1, jpk
            BaseMember(indexj,:,:,indexi)=e1t*e2t
        end do
        BaseMember(:,:,:,indexi)=BaseMember(:,:,:,indexi)*e3t_0*SeikMask
    end do
    
    MaxVarVec=MaxVarVec/BaseMember
    
    call CutCellsTracer(BaseMember)
    
    totalsum=0.0d0
    partialsum=sum(BaseMember(:,:,:,1))
    
    BaseMember=BaseMember*SeikTrcMask

    call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
    if (lwp) write (*,*) "total sum=", totalsum
    
    ModelErrorDiag1=ModelErrorDiag1/totalsum*reshape(BaseMember,(/SpaceDim/))

!call trcwriSeik("12345678901234567",-1,"RESTARTS/",BaseMember)
if (.false.) then !true to write a map of processes and coordinates
    BaseMember =Huge(BaseMember(1,1,1,1))
    BaseMember(:,:,:,1)=dble(myrank)
    do indexi=1, jpi
        BaseMember(:,:,indexi,2)=dble(indexi)
    end do
    do indexi=1, jpj
        BaseMember(:,indexi,:,3)=dble(indexi)
    end do
    call trcwriSeik("r="//trim(ctrcnm(1))//".x="//trim(ctrcnm(2))//".y="//trim(ctrcnm(3)),-1,"RESTARTS/",BaseMember)
end if

    BaseMember=Huge(BaseMember(1,1,1,1))
    
    ObsBaseMember=e1t*e2t*SeikMask(1,:,:)
    call CutCellsSurface(ObsBaseMember)
    
    totalsum=0.0d0
    partialsum=sum(ObsBaseMember)

    call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
    if (lwp) write (*,*) "total sum obs=", totalsum
    
    ObsErrorDiag1=ObsErrorDiag1/totalsum*reshape(ObsBaseMember,(/ObsSpaceDim/))
    
    ObsBaseMember=Huge(ObsBaseMember(1,1))
    
    if (UseDiffCov) then
        SeikUMask=(SeikMask(:,:,1:jpi-1)+SeikMask(:,:,2:jpi))/2
        SeikVMask=(SeikMask(:,1:jpj-1,:)+SeikMask(:,2:jpj,:))/2
        SeikWMask=(SeikMask(1:jpk-1,:,:)+SeikMask(2:jpk,:,:))/2
        
        SeikUTrcMask=(SeikTrcMask(:,:,1:jpi-1,:)+SeikTrcMask(:,:,2:jpi,:))/2
        SeikVTrcMask=(SeikTrcMask(:,1:jpj-1,:,:)+SeikTrcMask(:,2:jpj,:,:))/2
        SeikWTrcMask=(SeikTrcMask(1:jpk-1,:,:,:)+SeikTrcMask(2:jpk,:,:,:))/2
        
        do indexi=1, jptra
            do indexj=1, jpk
                UDiffBaseMember(indexj,:,:,indexi)=e1u(:,1:jpi-1)*e2u(:,1:jpi-1)
            end do
            UDiffBaseMember(:,:,:,indexi)=UDiffBaseMember(:,:,:,indexi)*e3u_0(:,:,1:jpi-1)*SeikUMask
        end do
        !if (CutLeft) UDiffBaseMember(:,:,1,:)=0.0d0
        if (CutRight) UDiffBaseMember(:,:,jpi-1,:)=0.0d0
        if (CutTop) UDiffBaseMember(:,jpj,:,:)=0.0d0
        if (CutBottom) UDiffBaseMember(:,1,:,:)=0.0d0
        
        totalsum=0.0d0
        partialsum=sum(UDiffBaseMember(:,:,:,1))
        UDiffBaseMember=UDiffBaseMember*SeikUTrcMask

        call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
        if (lwp) write (*,*) "total sum=", totalsum
        
        UDiffModelErrorDiag1=UDiffModelErrorDiag1/totalsum*reshape(UDiffBaseMember,(/UDiffSpaceDim/))!*16*10**12 !the med sea is about 4*10^6m long. I will use total sum of the obserror
        do indexi=1, jptra
            do indexj=1, jpk
                VDiffBaseMember(indexj,:,:,indexi)=e1v(1:jpj-1,:)*e2v(1:jpj-1,:)
            end do
            VDiffBaseMember(:,:,:,indexi)=VDiffBaseMember(:,:,:,indexi)*e3v_0(:,1:jpi-1,:)*SeikVMask
        end do
        if (CutLeft) VDiffBaseMember(:,:,1,:)=0.0d0
        if (CutRight) VDiffBaseMember(:,:,jpi,:)=0.0d0
        if (CutTop) VDiffBaseMember(:,jpj-1,:,:)=0.0d0
        !if (CutBottom) VDiffBaseMember(:,1,:,:)=0.0d0
        
        totalsum=0.0d0
        partialsum=sum(VDiffBaseMember(:,:,:,1))
        VDiffBaseMember=VDiffBaseMember*SeikVTrcMask
        
        call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
        if (lwp) write (*,*) "total sum=", totalsum
        
        VDiffModelErrorDiag1=VDiffModelErrorDiag1/totalsum*reshape(VDiffBaseMember,(/VDiffSpaceDim/))!*16*10**12 !the med sea is about 4*10^6m long
        
        do indexi=1, jptra
            do indexj=1, jpk-1
                WDiffBaseMember(indexj,:,:,indexi)=e1t*e2t
            end do
            WDiffBaseMember(:,:,:,indexi)=WDiffBaseMember(:,:,:,indexi)*e3w_0(2:jpk,:,:)*SeikWMask
        end do
        if (CutLeft) WDiffBaseMember(:,:,1,:)=0.0d0
        if (CutRight) WDiffBaseMember(:,:,jpi,:)=0.0d0
        if (CutTop) WDiffBaseMember(:,jpj,:,:)=0.0d0
        if (CutBottom) WDiffBaseMember(:,1,:,:)=0.0d0
                
        totalsum=0.0d0
        partialsum=sum(WDiffBaseMember(:,:,:,1))
        WDiffBaseMember=WDiffBaseMember*SeikWTrcMask
        
        call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
        if (lwp) write (*,*) "total sum=", totalsum
        
        WDiffModelErrorDiag1=WDiffModelErrorDiag1/totalsum*reshape(WDiffBaseMember,(/WDiffSpaceDim/))*25*10**6 !the med sea is about 5*10^3m deep. 
        !in this case it is convenient to set a different rescaling for horizontal and vertical covariance. probably it is advisable to keep WDiffModelErrorDiag1 small.
        
        UDiffObsBaseMember=e1u(:,1:jpi-1)*e2u(:,1:jpi-1)*SeikUMask(1,:,:)
        !if (CutLeft) UDiffObsBaseMember(:,1)=0.0d0
        if (CutRight) UDiffObsBaseMember(:,jpi-1)=0.0d0
        if (CutTop) UDiffObsBaseMember(jpj,:)=0.0d0
        if (CutBottom) UDiffObsBaseMember(1,:)=0.0d0
        
        totalsum=0.0d0
        partialsum=sum(UDiffObsBaseMember)

        call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
        if (lwp) write (*,*) "total sum obs=", totalsum
        
        UDiffObsErrorDiag1=UDiffObsErrorDiag1*reshape(UDiffObsBaseMember,(/UDiffObsSpaceDim/)) !/totalsum*16*10**12 !the med sea is about 4*10^6m long
        UDiffModelErrorDiag1=UDiffModelErrorDiag1*totalsum
        
        VDiffObsBaseMember=e1v(1:jpj-1,:)*e2v(1:jpj-1,:)*SeikVMask(1,:,:)
        if (CutLeft) VDiffObsBaseMember(:,1)=0.0d0
        if (CutRight) VDiffObsBaseMember(:,jpi)=0.0d0
        if (CutTop) VDiffObsBaseMember(jpj-1,:)=0.0d0
        !if (CutBottom) VDiffObsBaseMember(1,:)=0.0d0
        
        totalsum=0.0d0
        partialsum=sum(VDiffObsBaseMember)

        call mpi_allreduce(partialsum, totalsum, 1, MPI_real8, MPI_sum, LocalComm, ierr)
        if (lwp) write (*,*) "total sum obs=", totalsum
        
        VDiffObsErrorDiag1=VDiffObsErrorDiag1*reshape(VDiffObsBaseMember,(/VDiffObsSpaceDim/)) !/totalsum*16*10**12 !the med sea is about 4*10^6m long
        VDiffModelErrorDiag1=VDiffModelErrorDiag1*totalsum
        
        do indexi=1, UDiffSpaceDim
            if (.not.(UDiffModelErrorDiag1(indexi)==UDiffModelErrorDiag1(indexi))) then
                write(*,*) "myrank=" ,myrank, " indexi=", indexi
                error stop "NAN UDiffModelErrorDiag1"
            end if
        end do
        
        do indexi=1, VDiffSpaceDim
            if (.not.(VDiffModelErrorDiag1(indexi)==VDiffModelErrorDiag1(indexi))) then
                write(*,*) "myrank=" ,myrank, " indexi=", indexi
                error stop "NAN VDiffModelErrorDiag1"
            end if
        end do
        
        do indexi=1, WDiffSpaceDim
            if (.not.(WDiffModelErrorDiag1(indexi)==WDiffModelErrorDiag1(indexi))) then
                write(*,*) "myrank=" ,myrank, " indexi=", indexi
                error stop "NAN WDiffModelErrorDiag1"
            end if
        end do
        
        do indexi=1, UDiffObsSpaceDim
            if (.not.(UDiffObsErrorDiag1(indexi)==UDiffObsErrorDiag1(indexi))) then
                write(*,*) "myrank=" ,myrank, " indexi=", indexi
                error stop "NAN UDiffObsErrorDiag1"
            end if
        end do
        
        do indexi=1, VDiffObsSpaceDim
            if (.not.(VDiffObsErrorDiag1(indexi)==VDiffObsErrorDiag1(indexi))) then
                write(*,*) "myrank=" ,myrank, " indexi=", indexi
                error stop "NAN VDiffObsErrorDiag1"
            end if
        end do
        
        
    end if

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
