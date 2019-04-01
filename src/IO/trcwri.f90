      SUBROUTINE trcwri(datestring)


      USE myalloc
      USE IO_mem
      USE calendar
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module


      IMPLICIT NONE
      CHARACTER(LEN=17), INTENT(IN) :: datestring



!----------------------------------------------------------------------
! local declarations
! ==================
      double precision ::  Miss_val =1.e20
      INTEGER jk,jj,ji,jn
      double precision julian


      CHARACTER(LEN=37) filename

      CHARACTER(LEN=3) varname

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER ind1, i_contribution, j_contribution


       filename = 'RST.20111231-15:30:00.N1p.nc'

       julian=datestring2sec(datestring)

       if(lwp)write(*,*) 'trcwri ------------  myrank =',myrank,' datestring = ',  datestring

       trcwriparttime = MPI_WTIME() ! cronometer-start

      call mppsync()

      buf     = Miss_val
      bufftrn = Miss_val
      if (lwp) tottrn = Miss_val

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
                                    buf(jk,jj,ji) = trn(jk,jj,ji, jn)
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
                 bufftrn(ind1)= trn(jk,jj,ji, jn)
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

            varname=ctrcnm(jn)
            filename = 'RESTARTS/RST.'//datestring//'.'//varname//'.nc'

        CALL write_restart(filename,varname,julian, deflate_rst, deflate_level_rst)

        endif ! if myrank = 0
      END DO ! DO jn=1,jptra


       trcwriparttime = MPI_WTIME() - trcwriparttime
       trcwritottime = trcwritottime + trcwriparttime



      END SUBROUTINE trcwri
