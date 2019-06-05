!Sarebbe meglio riscrivere tutto qui... Ora va un po meglio, ma c'e' ancora da fare...

SUBROUTINE trcwriSeik(TimeString, BaseIndex, Directory, Tracer)
    use StringSEIK
    use myalloc
    
    implicit none
    CHARACTER(LEN=17), INTENT(IN) :: TimeString
    integer, intent(in) :: BaseIndex
    CHARACTER(len=*), intent(in) :: Directory
    double precision, dimension(jpk,jpj,jpi,jptra), intent(in) :: Tracer
    
    !Directory='REDUCED_BASE/DIMENSION_010/'
    !FileName = 'BASE001.20111231-15:30:00.N1p.nc'
    
    if (BaseIndex==-1) then
        call trcwriSeikWrapped(Directory//'RST', TimeString, Tracer)
    else
        call trcwriSeikWrapped(Directory//'BASE'//int2str(BaseIndex,3) , TimeString, Tracer)
    end if
end SUBROUTINE

subroutine trcwriSeikWrapped(PathAndFile, TimeString, Tracer)

      USE myalloc
      USE IO_mem
      USE calendar
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module


      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: PathAndFile
      CHARACTER(LEN=17), INTENT(IN) :: TimeString
      double precision, dimension(jpk,jpj,jpi,jptra), intent(in) :: Tracer
      
      


!----------------------------------------------------------------------
! local declarations
! ==================
      double precision ::  Miss_val =1.d20
      INTEGER jk,jj,ji,jn
      double precision julian

      CHARACTER(LEN=3) varname

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER ind1, i_contribution, j_contribution

        Miss_val =1.d20

       if (MyRank==0) write(*,*) 'trcwriSEIK ',PathAndFile

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
                                    buf(jk,jj,ji) = Tracer(jk,jj,ji, jn)
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
                 bufftrn(ind1)= Tracer(jk,jj,ji, jn)
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
            
            CALL write_restartSeik(PathAndFile//'.'//TimeString//'.'//varname//'.nc',varname,TimeString, deflate_rst, deflate_level_rst)

        endif ! if myrank = 0
      END DO ! DO jn=1,jptra


       trcwriparttime = MPI_WTIME() - trcwriparttime
       trcwritottime = trcwritottime + trcwriparttime



END SUBROUTINE 


SUBROUTINE Save2DSeik(datefrom,dateTo,PathAndFile, Tracer)
!     ******************
      !USE calendar
      USE myalloc
      USE IO_mem
      !USE FN_mem
      !USE TIME_MANAGER
      use mpi
      !USE ogstm_mpi_module

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) ::  dateFrom, dateTo
      character(len=*), intent(in) :: PathAndFile
      double precision, dimension(jpj, jpi), intent(in) :: Tracer
      
      
      INTEGER jj,ji
      INTEGER ind, i_contribution

      CHARACTER(LEN=20)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      double precision ::  Miss_val 

        Miss_val =1.e20
        if (myrank == 0) then
             tottrnIO2d = Miss_val
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

              tottrnIO2d (totjstart:totjend,totistart:totiend) = Tracer( reljstart:reljend,relistart:reliend)


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
            do ji =1 , jpi
                i_contribution = jpj * (ji-1)
                    do jj =1 , jpj
                    ind            = jj + i_contribution
                    buffDIA2d (ind)= Tracer(jj,ji)
                enddo
            enddo

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
              var        =  "MIS"

            d2f2d = REAL(tottrnIO2d(:,:),4)
            CALL WRITE_AVE_2dSeik(PathAndFile,trim(var),datefrom,dateTo, d2f2d)
 
      end if ! if(myrank == 0)
end subroutine















SUBROUTINE write_restartSeik(fileNetCDF,VAR, TimeString, deflate, deflate_level)
       USE netcdf
       USE myalloc

       IMPLICIT NONE
       CHARACTER*(*),intent(in) :: fileNetCDF
       CHARACTER(*),intent(in) ::  VAR
       CHARACTER(len=17), intent(in) :: TimeString
       integer, intent(in) :: deflate, deflate_level

       ! local
       integer :: istart, iend
       integer :: s, nc, counter
       integer :: timid, depid, yid, xid, xaid, yaid, zaid
       integer :: idB, idN, idLon, idLat, idLev, idTim
       integer shuffle
       
       integer indexi
       
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

       do indexi=1, jpjglo
            copy_inSeik(:,indexi,:)=transpose(tottrn(:,indexi,:))
       end do

       !call switch_index_double(tottrn,copy_inSeik,jpiglo,jpjglo,jpk)
       !call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idN,      copy_inSeik)
       call handle_err1(s,counter,fileNetCDF)
       
        s =nf90_close(nc)


END SUBROUTINE write_restartSeik

SUBROUTINE WRITE_AVE_2DSeik(fileNetCDF,VAR, datefrom, dateTo,M)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       character*(*),intent(in) :: fileNetCDF
       character(LEN=17),intent(in) :: datefrom, dateTo
       real,intent(in),dimension(jpjglo, jpiglo) :: M

       character(LEN=*) :: VAR
       integer :: istart,iend

       integer :: s, nc, counter
       integer :: timid, yid, xid
       integer :: idvartime,idphit,idlamt,idVAR
       real :: lat_actual_range(2), lon_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)



        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

       s = nf90_def_var(nc,trim(VAR) ,        nf90_float, (/xid,yid,timid/),  idVAR)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR, 'long_name'    ,VAR)
        s = nf90_put_att(nc,idVAR, 'missing_value' ,1.e+20)

        s =nf90_enddef(nc)

        counter=0
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,  transpose(M) )                    
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

END SUBROUTINE WRITE_AVE_2DSeik
