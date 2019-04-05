      SUBROUTINE diadump(datemean,datefrom,dateTo,FREQ_GROUP)
!     ******************
      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER jk,jj,ji, jn, jn_high
      INTEGER ind, i_contribution, j_contribution
      CHARACTER(10) newfile
      CHARACTER(LEN=42) forcing_file
      CHARACTER(LEN=60) bkpname
      CHARACTER(LEN=11) DIR
      logical IsBackup
      integer ave_counter


      CHARACTER(LEN=56) dia_file_nc
      CHARACTER(LEN=20)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      double precision ::  Miss_val =1.e20



     ! call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
! ----------------------------------------
      bkpname  = 'ave.20111231-15:30:00.N1p.nc.bkp'
      if (IsBackup) then
         forcing_file   = 'AVE_PHYS/ave.'//datemean//'.phys.nc.bkp'
      else
         forcing_file   = 'AVE_PHYS/ave.'//datemean//'.phys.nc'
      endif


      SELECT CASE (FREQ_GROUP)
        CASE (1) 
       ave_counter=ave_counter_1 
       DIR='AVE_FREQ_1/'
        CASE (2) 
       ave_counter=ave_counter_2 
       DIR='AVE_FREQ_2/'
      END SELECT

      if (myrank == 0) then !previously it was if lwp, but probably is wrong
          totsnIO   = Miss_val
          tottnIO   = Miss_val
          totunIO   = Miss_val
          totvnIO   = Miss_val
          totwnIO   = Miss_val
          totavtIO  = Miss_val
          tote3tIO  = Miss_val
          totvatmIO = Miss_val
          totempIO  = Miss_val
          totqsrIO  = Miss_val
      endif

!      PHYSICS FIRST!!
      if ( freq_ave_phys.eq.FREQ_GROUP) then
      ! *************** START COLLECTING DATA *****************************
       if (myrank == 0) then                    ! IF LABEL 1


! ******* myrank 0 sets indexes of tot matrix where to place its own part

          iPd    = nldi
          iPe    = nlei
          jPd    = nldj
          jPe    = nlej
          istart = nimpp
          jstart = njmpp
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

          totsnIO  (:, totjstart:totjend,totistart:totiend) = snIO    (:, reljstart:reljend,relistart:reliend)
          tottnIO  (:, totjstart:totjend,totistart:totiend) = tnIO    (:, reljstart:reljend,relistart:reliend)
          totvatmIO(totjstart:totjend,totistart:totiend) = vatmIO  (reljstart:reljend,relistart:reliend)
          totempIO (totjstart:totjend,totistart:totiend) = empIO   (reljstart:reljend,relistart:reliend)
          totqsrIO (totjstart:totjend,totistart:totiend) = qsrIO   (reljstart:reljend,relistart:reliend)
          totunIO  (:, totjstart:totjend,totistart:totiend) = unIO    (:, reljstart:reljend,relistart:reliend)
          totvnIO  (:, totjstart:totjend,totistart:totiend) = vnIO    (:, reljstart:reljend,relistart:reliend)
          totwnIO  (:, totjstart:totjend,totistart:totiend) = wnIO    (:, reljstart:reljend,relistart:reliend)
          totavtIO (:, totjstart:totjend,totistart:totiend) = avtIO   (:, reljstart:reljend,relistart:reliend)
          tote3tIO (:, totjstart:totjend,totistart:totiend) = e3tIO   (:, reljstart:reljend,relistart:reliend)

           do idrank = 1,CommSize-1
 ! **************  myrank 0 is receiving from the others their buffer  ****
               call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 1,LocalComm, status, ierr) !* first info to know where idrank is working
               call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 2,LocalComm, status, ierr)
               call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 3,LocalComm, status, ierr)
               call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 4,LocalComm, status, ierr)
               call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 5,LocalComm, status, ierr)
               call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 6,LocalComm, status, ierr)
               call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 7,LocalComm, status, ierr)
               call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 8,LocalComm, status, ierr)
       call MPI_RECV(buffsn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 11,LocalComm, status, ierr)
       call MPI_RECV(bufftn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 12,LocalComm, status, ierr)
       call MPI_RECV(buffvatm,jpi_rec*jpj_rec              ,mpi_real8,idrank, 13,LocalComm, status, ierr)
       call MPI_RECV(buffemp ,jpi_rec*jpj_rec              ,mpi_real8,idrank, 14,LocalComm, status, ierr)
       call MPI_RECV(buffqsr ,jpi_rec*jpj_rec              ,mpi_real8,idrank, 15,LocalComm, status, ierr)
       call MPI_RECV(buffun  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 16,LocalComm, status, ierr)
       call MPI_RECV(buffvn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 17,LocalComm, status, ierr)
       call MPI_RECV(buffwn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 18,LocalComm, status, ierr)
       call MPI_RECV(buffavt ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 19,LocalComm, status, ierr)
       call MPI_RECV(buffe3t ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 19,LocalComm, status, ierr)
 ! ******* myrank 0 sets indexes of tot matrix where to place buffers of idrank
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
               do ji =totistart,totiend ! 3d vars
                    i_contribution = jpk*jpj_rec*(ji-1 - totistart+ relistart )
                     do jj =totjstart,totjend
                         j_contribution = jpk*(jj-totjstart+ reljstart-1)
                         do jk =1 , jpk
                            ind = jk + j_contribution + i_contribution
                            totsnIO (jk,jj,ji)= buffsn (ind)
                            tottnIO (jk,jj,ji)= bufftn (ind)
                            totunIO (jk,jj,ji)= buffun (ind)
                            totvnIO (jk,jj,ji)= buffvn (ind)
                            totwnIO (jk,jj,ji)= buffwn (ind)
                            totavtIO(jk,jj,ji)= buffavt(ind)
                            tote3tIO(jk,jj,ji)= buffe3t(ind)
                         enddo
                      enddo
                   enddo

             do ji =totistart,totiend  ! and 2d vars
                i_contribution = jpj_rec*(ji-1 -totistart+ relistart )
                do jj =totjstart,totjend
                   ind = (jj - totjstart+ reljstart) + i_contribution
                   totvatmIO (jj,ji)= buffvatm(ind)
                   totempIO  (jj,ji)= buffemp (ind)
                   totqsrIO  (jj,ji)= buffqsr (ind)
                enddo
             enddo
           enddo !idrank = 1, size-1
       else  ! IF LABEL 1,  if(myrank == 0)
            do ji =1 , jpi
             i_contribution= jpk*jpj * (ji - 1 )
             do jj =1 , jpj
              j_contribution=jpk*(jj-1)
              do jk =1 , jpk
                    ind         =  jk + j_contribution + i_contribution
                    buffsn (ind)= snIO (jk,jj,ji)
                    bufftn (ind)= tnIO (jk,jj,ji)
                    buffun (ind)= unIO (jk,jj,ji)
                    buffvn (ind)= vnIO (jk,jj,ji)
                    buffwn (ind)= wnIO (jk,jj,ji)
                    buffavt(ind)= avtIO(jk,jj,ji)
                    buffe3t(ind)= e3tIO(jk,jj,ji)
               enddo
              enddo
             enddo
             do ji =1 , jpi
              i_contribution=jpj * (ji-1)
              do jj =1 , jpj
                ind           = jj + i_contribution
                buffvatm (ind)= vatmIO(jj,ji)
                buffemp  (ind)= empIO (jj,ji)
                buffqsr  (ind)= qsrIO (jj,ji)
               enddo
              enddo
               call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, LocalComm,ierr)
               call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, LocalComm,ierr)
               call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, LocalComm,ierr)
               call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, LocalComm,ierr)
               call MPI_SEND(nlei , 1,mpi_integer, 0, 5, LocalComm,ierr)
               call MPI_SEND(nlej , 1,mpi_integer, 0, 6, LocalComm,ierr)
               call MPI_SEND(nldi , 1,mpi_integer, 0, 7, LocalComm,ierr)
               call MPI_SEND(nldj , 1,mpi_integer, 0, 8, LocalComm,ierr)
            call MPI_SEND(buffsn  , jpk*jpj*jpi  ,mpi_real8, 0, 11, LocalComm,ierr)
            call MPI_SEND(bufftn  , jpk*jpj*jpi  ,mpi_real8, 0, 12, LocalComm,ierr)
            call MPI_SEND(buffvatm, jpi*jpj      ,mpi_real8, 0, 13, LocalComm,ierr)
            call MPI_SEND(buffemp , jpi*jpj      ,mpi_real8, 0, 14, LocalComm,ierr)
            call MPI_SEND(buffqsr , jpi*jpj      ,mpi_real8, 0, 15, LocalComm,ierr)
            call MPI_SEND(buffun  , jpk*jpj*jpi  ,mpi_real8, 0, 16, LocalComm,ierr)
            call MPI_SEND(buffvn  , jpk*jpj*jpi  ,mpi_real8, 0, 17, LocalComm,ierr)
            call MPI_SEND(buffwn  , jpk*jpj*jpi  ,mpi_real8, 0, 18, LocalComm,ierr)
            call MPI_SEND(buffavt , jpk*jpj*jpi  ,mpi_real8, 0, 19, LocalComm,ierr)
            call MPI_SEND(buffe3t , jpk*jpj*jpi  ,mpi_real8, 0, 19, LocalComm,ierr)




       endif ! IF LABEL 1, if(myrank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************

      if(myrank == 0) then ! IF LABEL 4,
         if (IsBackup) then

            call PhysDump_bkp(forcing_file, datefrom, dateTo,ave_counter)
          else
            call PhysDump(forcing_file, datefrom, dateTo)
         endif
      endif

      endif !if ( freq_ave_phys.eq.FREQ_GROUP)



      jn_high = 0
      if (myrank == 0 ) tottrnIO2d = Miss_val !previously it was if lwp, but probably is wrong
! ******************  DIAGNOSTIC OUTPUT   2D *******************
      DO jn = 1, JPTRA_DIA_2D

           if (.not.is_time_to_save(jn,FREQ_GROUP,2)) CYCLE
           if (FREQ_GROUP.eq.1) jn_high = jn_high+1
       if (myrank == 0) then
!   ! ******* myrank 0 sets indexes of tot matrix where to place its own part

             iPd    = nldi
             iPe    = nlei
             jPd    = nldj
             jPe    = nlej
             istart = nimpp
             jstart = njmpp
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

              if (FREQ_GROUP.eq.1) then
              tottrnIO2d (totjstart:totjend,totistart:totiend) = tra_DIA_2d_IO_HIGH(jn_high, reljstart:reljend,relistart:reliend)
              else
              tottrnIO2d (totjstart:totjend,totistart:totiend) = tra_DIA_2d_IO(jn,  reljstart:reljend,relistart:reliend)
              endif

              do idrank = 1,CommSize-1
 ! **************  myrank 0 is receiving from the others their buffer  ****
                 call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 32,LocalComm, status, ierr) !* first info to know where idrank is working
                 call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 33,LocalComm, status, ierr)
                 call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 34,LocalComm, status, ierr)
                 call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 35,LocalComm, status, ierr)
                 call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 36,LocalComm, status, ierr)
                 call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 37,LocalComm, status, ierr)
                 call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 38,LocalComm, status, ierr)
                 call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 39,LocalComm, status, ierr)
                 call MPI_RECV(buffDIA2d,jpi_rec*jpj_rec      ,mpi_real8,idrank, 40,LocalComm, status, ierr)
 ! ******* myrank 0 sets indexes of tot matrix where to place buffers of idrank
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
               do ji =totistart,totiend ! only 2d vars
                 i_contribution = jpj_rec*(ji-totistart+ relistart -1)
                 do jj =totjstart,totjend
                   ind = jj-totjstart+ reljstart + i_contribution
                   tottrnIO2d (jj,ji)= buffDIA2d(ind)
                  enddo
                 enddo
              enddo !idrank = 1, size-1
       ELSE ! ranks 1 --> size-1
       if (FREQ_GROUP.eq.2) then
             do ji =1 , jpi
              i_contribution = jpj * (ji-1)
              do jj =1 , jpj
                ind            = jj + i_contribution
                buffDIA2d (ind)= tra_DIA_2d_IO(jn,jj,ji)
               enddo
              enddo
       else
             do ji =1 , jpi
              i_contribution = jpj * (ji-1)
              do jj =1 , jpj
                ind           = jj + i_contribution
                buffDIA2d (ind)= tra_DIA_2d_IO_high(jn_high,jj,ji)
               enddo
              enddo
       endif
               call MPI_SEND(jpi  , 1,mpi_integer, 0, 32, LocalComm,ierr)
               call MPI_SEND(jpj  , 1,mpi_integer, 0, 33, LocalComm,ierr)
               call MPI_SEND(nimpp, 1,mpi_integer, 0, 34, LocalComm,ierr)
               call MPI_SEND(njmpp, 1,mpi_integer, 0, 35, LocalComm,ierr)
               call MPI_SEND(nlei , 1,mpi_integer, 0, 36, LocalComm,ierr)
               call MPI_SEND(nlej , 1,mpi_integer, 0, 37, LocalComm,ierr)
               call MPI_SEND(nldi , 1,mpi_integer, 0, 38, LocalComm,ierr)
               call MPI_SEND(nldj , 1,mpi_integer, 0, 39, LocalComm,ierr)
              call MPI_SEND(buffDIA2d, jpi*jpj   ,mpi_real8, 0, 40, LocalComm,ierr)
       ENDIF

      if (myrank == 0) then
              var        =  dianm_2d(jn)
              bkpname     = DIR//'ave.'//datemean//'.'//trim(var)//'.nc.bkp'
              dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var)//'.nc'

              if (IsBackup) then
                 !write(*,*) "trcdia ave_counter --> ", bkpname, ave_counter
                 CALL WRITE_AVE_2d_BKP(bkpname,var,datefrom, dateTo,tottrnIO2d, ave_counter)
      
              else
                 d2f2d = REAL(tottrnIO2d(:,:),4)
                 CALL WRITE_AVE_2d(dia_file_nc,var,datefrom,dateTo, d2f2d)
      

              endif
      end if ! if(myrank == 0)

         if (.not.IsBackup) then
             if (FREQ_GROUP.eq.2) then
                tra_DIA_2d_IO(jn,:,:) = 0.
              else
                tra_DIA_2d_IO_HIGH(jn_high,:,:) = 0.
              endif
          endif

      ENDDO ! on jn


       jn_high = 0
! ! ******************  3D DIAGNOSTIC OUTPUT   *******************
       DO jn =1 , jptra_dia

           if (.not.is_time_to_save(jn,FREQ_GROUP,3)) CYCLE
           if (FREQ_GROUP.eq.1) jn_high = jn_high+1

           if (myrank == 0) then                    ! IF LABEL 1


! ! ******* myrank 0 sets indexes of tot matrix where to place its own part

             iPd    = nldi
             iPe    = nlei
             jPd    = nldj
             jPe    = nlej
             istart = nimpp
             jstart = njmpp
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

      if (FREQ_GROUP.eq.1) then
      tottrnIO (:,totjstart:totjend,totistart:totiend) = tra_DIA_IO_HIGH(jn_high,:,reljstart:reljend,relistart:reliend)
      else
      tottrnIO (:, totjstart:totjend,totistart:totiend) = tra_DIA_IO(jn,:, reljstart:reljend, relistart:reliend) ! diagnostic from reaction model
      endif
              do idrank = 1,CommSize-1
 ! **************  myrank 0 is receiving from the others their buffer  ****
                 call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 22,LocalComm, status, ierr) !* first info to know where idrank is working
                 call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 23,LocalComm, status, ierr)
                 call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 24,LocalComm, status, ierr)
                 call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 25,LocalComm, status, ierr)
                 call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 26,LocalComm, status, ierr)
                 call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 27,LocalComm, status, ierr)
                 call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 28,LocalComm, status, ierr)
                 call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 29,LocalComm, status, ierr)
                 call MPI_RECV(buffDIA  ,jpi_rec*jpj_rec*jpk  ,mpi_real8   ,idrank, 30,LocalComm, status, ierr)
 ! ******* myrank 0 sets indexes of tot matrix where to place buffers of idrank
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
                 do ji =totistart,totiend ! 3d vars
                    i_contribution = jpk*jpj_rec*(ji-1 -totistart+ relistart )
                    do jj =totjstart,totjend
                       j_contribution = jpk*(jj-totjstart+ reljstart-1)
                       do jk =1 , jpk
                            ind = jk + j_contribution + i_contribution
                            tottrnIO(jk,jj,ji)= buffDIA (ind)
                       enddo
                    enddo
                 enddo
              enddo !idrank = 1, size-1
           else  ! IF LABEL 1,  if(myrank == 0)
       if (FREQ_GROUP.eq.2) then
             do ji =1, jpi
               i_contribution= jpk*jpj * (ji - 1 )
               do jj =1 , jpj
                 j_contribution=jpk*(jj-1)
                 do jk =1 , jpk
                   ind         =  jk + j_contribution + i_contribution
                   buffDIA(ind) = tra_DIA_IO(jn, jk,jj,ji)
                 enddo
               enddo
             enddo
       else
             do ji =1, jpi
               i_contribution= jpk*jpj * (ji - 1 )
               do jj =1 , jpj
                 j_contribution=jpk*(jj-1)
                 do jk =1 , jpk
                   ind         =  jk + j_contribution + i_contribution
                   buffDIA(ind) = tra_DIA_IO_HIGH(jn_high, jk,jj,ji)
                 enddo
               enddo
             enddo
       endif
               call MPI_SEND(jpi  , 1,mpi_integer, 0, 22, LocalComm,ierr)
               call MPI_SEND(jpj  , 1,mpi_integer, 0, 23, LocalComm,ierr)
               call MPI_SEND(nimpp, 1,mpi_integer, 0, 24, LocalComm,ierr)
               call MPI_SEND(njmpp, 1,mpi_integer, 0, 25, LocalComm,ierr)
               call MPI_SEND(nlei , 1,mpi_integer, 0, 26, LocalComm,ierr)
               call MPI_SEND(nlej , 1,mpi_integer, 0, 27, LocalComm,ierr)
               call MPI_SEND(nldi , 1,mpi_integer, 0, 28, LocalComm,ierr)
               call MPI_SEND(nldj , 1,mpi_integer, 0, 29, LocalComm,ierr)
               call MPI_SEND(buffDIA  , jpk*jpj*jpi,mpi_real8, 0, 30, LocalComm,ierr)
       endif ! IF LABEL 1, if(myrank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************
      if (myrank == 0) then
              var        =  dianm(jn)

              bkpname     = DIR//'ave.'//datemean//'.'//trim(var)//'.nc.bkp'
              dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var)//'.nc'
              
              if (IsBackup) then
                 CALL WRITE_AVE_BKP(bkpname,var,datefrom, dateTo,tottrnIO,ave_counter,deflate_ave, deflate_level_ave)
              else
                 CALL WRITE_AVE(dia_file_nc,var,datefrom,dateTo, tottrnIO,deflate_ave, deflate_level_ave)
              endif


      end if ! IF LABEL 4  if(myrank == 0)
         if (.not.IsBackup) then
             if (FREQ_GROUP.eq.2) then
                tra_DIA_IO(jn,:,:,:) = 0.
              else
                tra_DIA_IO_HIGH(jn_high,:,:,:) = 0.
              endif
          endif
      enddo  ! loop in jn



      if ((.not.IsBackup).and.( freq_ave_phys.eq.FREQ_GROUP) ) then

          snIO     = 0.
          tnIO     = 0.
          vatmIO   = 0.
          empIO    = 0.
          qsrIO    = 0.
          unIO     = 0.
          vnIO     = 0.
          wnIO     = 0.
          avtIO    = 0.
          e3tIO    = 0
      endif


      CONTAINS

      LOGICAL FUNCTION IS_TIME_TO_SAVE(jn,FREQ_GROUP,ndims)
      IMPLICIT NONE

      integer jn, FREQ_GROUP,ndims

      IF (FREQ_GROUP.eq.2) then
         IS_TIME_TO_SAVE = .true.
      ELSE
        IS_TIME_TO_SAVE = .false.
        IF (ndims==3) then
            IF (diahf(jn).eq.1)   IS_TIME_TO_SAVE = .true.
        ELSE
           IF (diahf_2d(jn).eq.1) IS_TIME_TO_SAVE = .true.
        ENDIF
      ENDIF
      END FUNCTION IS_TIME_TO_SAVE

      end SUBROUTINE diadump
