!Sarebbe meglio riscrivere tutto qui...

SUBROUTINE trcwriSeik(TimeString, BaseIndex, Directory)


      USE myalloc
      USE IO_mem
      USE calendar
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module


      IMPLICIT NONE
      CHARACTER(LEN=17), INTENT(IN) :: TimeString
      integer, intent(in) :: BaseIndex
      CHARACTER(len=*), intent(in) :: Directory



!----------------------------------------------------------------------
! local declarations
! ==================
      double precision ::  Miss_val =1.d20
      INTEGER jk,jj,ji,jn
      double precision julian


      CHARACTER(LEN=32) FileName

      CHARACTER(LEN=3) varname

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER ind1, i_contribution, j_contribution

        Miss_val =1.d20
        !Directory='REDUCED_BASE/DIMENSION_010/'
        FileName = 'BASE001.20111231-15:30:00.N1p.nc'

       if(lwp)write(*,*) 'trcwriSEIK ------------  myrank =',myrank, ' BaseIndex=' , BaseIndex

       trcwriparttime = MPI_WTIME() ! cronometer-start

      call mppsync()

      buf     = Miss_val
      bufftrn = Miss_val
      if (myrank == 0) tottrn = Miss_val !previously it was lwp, but probably incorrect

       DO jn=1,jptra
        if(myrank == 0) then
           istart = nimpp
           jstart = njmpp
           iPd = nldi
           iPe = nlei
           jPd = nldj
           jPe = nlej
           irange    = iPe - iPd + 1
           jrange    = jPe - jPd + 1
           totistart = istart + iPd - 1
           totiend   = totistart + irange - 1
           totjstart = jstart + jPd - 1
           totjend   = totjstart + jrange - 1
           relistart = 1 + iPd - 1
           reliend   = relistart + irange - 1
           reljstart = 1 + jPd - 1
           reljend   = reljstart + jrange - 1


            do ji =1 , jpi
                  do jj =1 , jpj
                        do jk =1 , jpk
                              if (tmask(jk,jj,ji).eq.1) then
                                    buf(jk,jj,ji) = BaseMemeber(jk,jj,ji, jn)
                              endif
                        enddo
                  enddo
            enddo
           tottrn  (:,totjstart:totjend,totistart:totiend)= buf(:, reljstart:reljend,relistart:reliend)




      do idrank = 1,CommSize-1
         call MPI_RECV(jpi_rec , 1,                  mpi_integer, idrank,  1,LocalComm, status, ierr)
         call MPI_RECV(jpj_rec , 1,                  mpi_integer, idrank,  2,LocalComm, status, ierr)
         call MPI_RECV(istart  , 1,                  mpi_integer, idrank,  3,LocalComm, status, ierr)
         call MPI_RECV(jstart  , 1,                  mpi_integer, idrank,  4,LocalComm, status, ierr)
         call MPI_RECV(iPe     , 1,                  mpi_integer, idrank,  5,LocalComm, status, ierr)
         call MPI_RECV(jPe     , 1,                  mpi_integer, idrank,  6,LocalComm, status, ierr)
         call MPI_RECV(iPd     , 1,                  mpi_integer, idrank,  7,LocalComm, status, ierr)
         call MPI_RECV(jPd     , 1,                  mpi_integer, idrank,  8,LocalComm, status, ierr)
         call MPI_RECV(bufftrn,   jpi_rec*jpj_rec*jpk, mpi_real8, idrank, 11,LocalComm, status, ierr)



           irange    = iPe - iPd + 1
           jrange    = jPe - jPd + 1
           totistart = istart + iPd - 1
           totiend   = totistart + irange - 1
           totjstart = jstart + jPd - 1
           totjend   = totjstart + jrange - 1
           relistart = 1 + iPd - 1
           reliend   = relistart + irange - 1
           reljstart = 1 + jPd - 1
           reljend   = reljstart + jrange - 1

          do ji =totistart,totiend
            i_contribution   = jpk*jpj_rec*(ji-1-totistart+ relistart)
            do jj =totjstart,totjend
              j_contribution = jpk*(jj-1-totjstart+ reljstart)
              do jk =1, jpk
                  ind1 = jk + j_contribution + i_contribution
                  tottrn(jk,jj,ji)= bufftrn(ind1)
              enddo
            enddo
          enddo

          enddo ! do idrank = 1, size-1
       else !myrank != 0
           do ji =1 , jpi
            i_contribution= jpk*jpj * (ji - 1 )
            do jj =1 , jpj
             j_contribution=jpk*(jj-1)
             do jk =1 , jpk
              ind1 = jk + j_contribution + i_contribution
              if (tmask(jk,jj,ji).eq.1) then
                 bufftrn(ind1)= BaseMemeber(jk,jj,ji, jn)
              endif
             enddo
            enddo
           enddo
           call MPI_SEND(jpi      , 1         ,mpi_integer, 0,  1, LocalComm,ierr)
           call MPI_SEND(jpj      , 1         ,mpi_integer, 0,  2, LocalComm,ierr)
           call MPI_SEND(nimpp    , 1         ,mpi_integer, 0,  3, LocalComm,ierr)
           call MPI_SEND(njmpp    , 1         ,mpi_integer, 0,  4, LocalComm,ierr)
           call MPI_SEND(nlei     , 1         ,mpi_integer, 0,  5, LocalComm,ierr)
           call MPI_SEND(nlej     , 1         ,mpi_integer, 0,  6, LocalComm,ierr)
           call MPI_SEND(nldi     , 1         ,mpi_integer, 0,  7, LocalComm,ierr)
           call MPI_SEND(nldj     , 1         ,mpi_integer, 0,  8, LocalComm,ierr)
           call MPI_SEND(bufftrn  ,jpk*jpj*jpi,  mpi_real8, 0, 11, LocalComm,ierr)
       endif ! if myrank = 0


        if(myrank == 0) then

            varname=trim(ctrcnm(jn))
            write(FileName,'(A4,I3.3,A25)') 'BASE', BaseIndex, '.'//TimeString//'.'//varname//'.nc'
            
            CALL write_restartSeik(Directory//FileName,varname,TimeString, deflate_rst, deflate_level_rst)

        endif ! if myrank = 0
      END DO ! DO jn=1,jptra


       trcwriparttime = MPI_WTIME() - trcwriparttime
       trcwritottime = trcwritottime + trcwriparttime



END SUBROUTINE trcwriSeik

SUBROUTINE write_restartSeik(fileNetCDF,VAR, TimeString, deflate, deflate_level)
       USE netcdf
       USE myalloc

       IMPLICIT NONE
       CHARACTER*(*),intent(in) :: fileNetCDF
       CHARACTER(*),intent(in) ::  VAR
       CHARACTER(len=17) intent(in) :: TimeString
       integer, intent(in) :: deflate, deflate_level

       ! local
       integer :: istart, iend
       integer :: s, nc, counter
       integer :: timid, depid, yid, xid, xaid, yaid, zaid
       integer :: idB, idN, idLon, idLat, idLev, idTim
       integer shuffle
       
       !TimeString =fileNetCDF(14:30)
       shuffle       = 0

       ! Just to try without 'or'
       ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
       s = nf90_create(fileNetCDF, NF90_HDF5, nc)

      s = nf90_put_att(nc, nf90_global, 'TimeString'     , TimeString)
        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'   , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'   , jpjglo,  yid)
        s= nf90_def_dim(nc,'z'   , jpk   ,depid)
        s= nf90_def_dim(nc,'time', 1     ,timid)
        s= nf90_def_dim(nc,'x_a'      , 1,  xaid)
        s= nf90_def_dim(nc,'y_a'      , 1,  yaid)
        s= nf90_def_dim(nc,'z_a'      , 3  ,zaid)


       s = nf90_def_var(nc,'nav_lon', nf90_double,  (/xid,yid/), idLon)
       s = nf90_def_var(nc,'nav_lat', nf90_double,  (/xid,yid/), idLat)
       s = nf90_def_var(nc,'nav_lev', nf90_double,  (/depid/)  , idLev)
       !s = nf90_def_var(nc,'time'   , nf90_double,  (/timid/)  , idTim)
        s = nf90_def_var(nc,'TRN'//VAR, nf90_double, (/xid,yid,depid,timid/), idN)
        s = nf90_def_var_deflate(nc, idN, shuffle, deflate, deflate_level)
        call handle_err1(s,counter,fileNetCDF)
        !s= nf90_put_att(nc,idTim ,'Units', 'seconds since 1582-10-15 00:00:00')
      
        s = nf90_put_att(nc,idN   , 'missing_value',1.d20)
        s =nf90_enddef(nc)
        s = nf90_put_var(nc, idLon,  TRANSPOSE(totglamt))
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idLat,  TRANSPOSE(totgphit))
       call handle_err1(s,counter,fileNetCDF)
       

        s = nf90_put_var(nc, idLev,     gdept)
       call handle_err1(s,counter,fileNetCDF)

       

       call switch_index_double(tottrn,copy_inSeik,jpiglo,jpjglo,jpk)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idN,      copy_inSeik)
       call handle_err1(s,counter,fileNetCDF)
       
        s =nf90_close(nc)


END SUBROUTINE write_restartSeik
