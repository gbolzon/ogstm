      SUBROUTINE trcrstDA(datestring)
!---------------------------------------------------------------------
!
!                       ROUTINE trcrst
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ files for restart for passive tracer
!
!----------------------------------------------------------------------


       USE calendar
       USE myalloc
       USE myalloc_mpp
       USE TIME_MANAGER
       USE DA_mem

       IMPLICIT NONE
      CHARACTER(LEN=17), INTENT(IN) :: datestring
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn, jDA
      CHARACTER(LEN=37) filename
      CHARACTER(LEN=6)  varname
      CHARACTER(LEN=43) bkpname




      DO jn=1, jptra  ! global loop on tracers to read restart


      if (.not.isaDAvar(ctrcnm(jn))) CYCLE


         varname  = 'TRN'//ctrcnm(jn)
         filename = 'RESTARTS/RST.'//datestring//'.'//trim(ctrcnm(jn))//'.nc'
         if (lwp) write(*,*) 'reading ', filename
         CALL readnc_slice_floatDA(filename,varname, trn(:,:,1:AssimilationLevels,jn) )



! ********************   we put initial undef to 0
          trb(:,:,:,jn) = trn(:,:,:,jn) * tmask; !! ACHTUNG !!!
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask;
      ENDDO

      END SUBROUTINE trcrstDA


      SUBROUTINE readnc_slice_floatDA(fileNetCDF,varname, M)
      USE myalloc
      USE myalloc_mpp
      USE netcdf
      USE DA_mem, ONLY : AssimilationLevels
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,AssimilationLevels)
      real(4), ALLOCATABLE :: Mfloat(:,:,:)
      
      ALLOCATE(Mfloat(jpi,jpj,AssimilationLevels))

      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, AssimilationLevels, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,Mfloat,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ;  call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      M(:,:,:) = real(Mfloat(:,:,:),8)
      DEALLOCATE(Mfloat)

      END SUBROUTINE readnc_slice_floatDA



