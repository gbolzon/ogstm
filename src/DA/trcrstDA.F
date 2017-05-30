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
      CHARACTER(LEN=15), INTENT(IN) :: datestring
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn, jDA
      integer, parameter :: jpk_200=26
      CHARACTER(LEN=37) filename
      CHARACTER(LEN=6)  varname
      CHARACTER(LEN=43) bkpname




      DO jn=1, jptra  ! global loop on tracers to read restart


      if (.not.isaDAvar(ctrcnm(jn))) CYCLE


         varname  = 'TRN'//ctrcnm(jn)
         filename = 'RESTARTS/RST.'//datestring//'.'//trim(ctrcnm(jn))//'.nc'
         if (lwp) write(*,*) 'reading ', filename
         CALL readnc_slice_doubleDA(filename,varname, trn(:,:,1:jpk_200,jn) )



! ********************   we put initial undef to 0
          trb(:,:,:,jn) = trn(:,:,:,jn) * tmask; !! ACHTUNG !!!
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask;
      ENDDO

      END SUBROUTINE trcrstDA


      SUBROUTINE readnc_slice_doubleDA(fileNetCDF,varname, M)
      USE myalloc
      USE myalloc_mpp
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid
      integer, parameter :: jpk_200=26

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk_200)

      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk_200, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ;  call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      END SUBROUTINE readnc_slice_doubleDA



