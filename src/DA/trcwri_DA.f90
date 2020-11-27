SUBROUTINE trcwriDA(datestring)
        !---------------------------------------------------------------------
        !
        !                       ROUTINE trcwriDA
        !
        !                     ******************
        !  gcoidessa@inogs.it  developing
        !
        !  Purpose :
        !  ---------
        !     Standard output of DA RESTARTS
        !  SUBROUTINE:
        !       - trcwriDA
        !       - write_BeforeAss
                
    
   

        USE netcdf
        USE myalloc
        USE IO_mem
        USE calendar
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module
        USE DA_mem

        USE MPI_GATHER_INFO

        USE MATRIX_VARS
        USE NODES_MODULE
        USE DTYPE_PROCS_STRING_MODULE
        USE DA_VARS_module
        
        IMPLICIT NONE
        CHARACTER(LEN=17), INTENT(IN) :: datestring



!----------------------------------------------------------------------
! local declarations
! ==================
        double precision ::  Miss_val =1.e20
        INTEGER jk,jj,ji,jn,i
        INTEGER s, nc, counter
        integer timid, depid, yid, xid, idvar
        double precision julian


        CHARACTER(LEN=45) BeforeName
        CHARACTER(LEN=43) BeforeNameShort

        CHARACTER(LEN=3) varname

        INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
        INTEGER irange, jrange
        INTEGER totistart, totiend, relistart, reliend
        INTEGER totjstart, totjend, reljstart, reljend
        INTEGER ind1, i_contribution, j_contribution
        INTEGER SysErr, system
        INTEGER :: jv,n_dumping_cycles,writing_rank,counter_var_DA,ivar,jn_da,ind_col
        CHARACTER(LEN=20)  var_to_store
        DOUBLE PRECISION :: packing_trwrida_sum_time,start_time_trcwrida,packing_trcwrida_init_time,gatherv_trcwrida_init_time,gatherv_trcwrida_fin_time,packing_trcwrida_fin_time,gatherv_trcwrida_delta_time,packing_trcwrida_delta_time,gatherv_trcwrida_sum_time,unpacking_trcwrida_rank_init_time,unpacking_trcwrida_rank_fin_time,writing_trcwrida_rank_fin_time,writing_trcwrida_rank_delta_time,unpacking_trcwrida_rank_delta_time,writing_trcwrida_rank_sum_time,unpacking_trcwrida_rank_sum_time,writing_trcwrida_rank_init_time,packing_trcwrida_sum_time
        

        julian=datestring2sec(datestring)

        if(WRITING_rank_wr)write(*,*) 'trcwri DA ------------  myrank =', myrank,' datestring = ', datestring 
        start_time_trcwrida = MPI_Wtime()       

        trcwriparttime = MPI_WTIME() ! cronometer-start
        call mppsync()
   
        buf     = Miss_val
        bufftrn = Miss_val
        if (WRITING_RANK_WR) tottrn = Miss_val
        n_dumping_cycles = matrix_DA_row
        counter_var_DA = 1
        writing_trcwrida_rank_sum_time = 0
        unpacking_trcwrida_rank_sum_time = 0
        packing_trwrida_sum_time = 0
        gatherv_trcwrida_fin_time = 0
        !all ranks
 
        DA_LOOP: DO jv = 1, n_dumping_cycles

<<<<<<< HEAD
                DO ivar = 1 , nodes!number of variables for each round corresponds to the number of nodes
                        packing_trcwrida_init_time = MPI_Wtime()
=======
                DO ivar = 1 , nodes*num_of_wr_procs_perNODE!number of variables for each round corresponds to the number of nodes

>>>>>>> 2a6e2ff45c1648244bcf74b65b8928f989d0c3c4
                        writing_rank = writing_procs(ivar)

                        IF (COUNTER_VAR_DA + PX_DA > num_DA_vars )then
                                EXIT
                        ELSE
                                jn_da = DA_table(COUNTER_VAR_DA + PX_DA)
                                !PACKING
                                do ji =1 , jpi
                                        i_contribution= jpk*jpj * (ji - 1 )
                                        do jj =1 , jpj
                                                j_contribution=jpk*(jj-1)
                                                do jk =1 , jpk
                                                        ind1 = jk + j_contribution + i_contribution
                                                        if (tmask(jk,jj,ji).eq.1) then
                                                                bufftrn(ind1)= trn(jk,jj,ji, jn_da)
                                                        endif
                                                enddo
                                        enddo
                                enddo

                                counter_var_DA = counter_var_DA + 1
                                gatherv_trcwrida_init_time = MPI_Wtime()
                                CALL MPI_GATHERV(bufftrn, sendcount, MPI_DOUBLE_PRECISION, bufftrn_TOT, jprcv_count, jpdispl_count, MPI_DOUBLE_PRECISION, writing_rank, MPI_COMM_WORLD, IERR)
                                gatherv_trcwrida_fin_time = gatherv_trcwrida_fin_time + MPI_Wtime()
                        END IF

                END DO
                packing_trcwrida_fin_time = MPI_Wtime()

                gatherv_trcwrida_delta_time = gatherv_trcwrida_fin_time - gatherv_trcwrida_init_time
                packing_trcwrida_delta_time = packing_trcwrida_fin_time - packing_trcwrida_init_time
                packing_trwrida_sum_time = packing_trwrida_sum_time + packing_trcwrida_delta_time

                CALL MPI_Reduce( gatherv_trcwrida_delta_time, gatherv_trcwrida_sum_time, 1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,IERROR)
                CALL MPI_Reduce( packing_trcwrida_delta_time, packing_trwrida_sum_time,2,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,IERROR)

                if(myrank==0)then
                        write(*,*)' gatherv_sum_time_trcwrida is ',gatherv_trcwrida_sum_time
                        write(*,*)'packing_sum_time_trcwrida is ',packing_trcwrida_sum_time
                end if



        !WRITING RANKS

                if(WRITING_RANK_WR) then
                                
                        do i=1, nodes*num_of_wr_procs_perNODE
                                if(myrank == writing_procs(i))then
                                        ind_col=i
                                        exit
                                end if
                        end do 
                
!                        ind_col = (myrank / n_ranks_per_node)+1

                        var_to_store = matrix_DA(jv,ind_col)%var_name
                        unpacking_trcwrida_rank_init_time = MPI_WTIME()
                        IF (var_to_store == "novars_input")then
                                EXIT
                        ELSE
                                DO idrank = 0,mpi_glcomm_size-1

                                        ! ******* WRITING RANK sets indexes of tot matrix where to place buffers of idrank
                                        irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange - 1
                                        totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange - 1
                                        relistart = 1 + iPd_a(idrank+1) - 1
                                        reliend   = relistart + irange - 1
                                        reljstart = 1 + jPd_a(idrank+1) - 1
                                        reljend   = reljstart + jrange - 1

                                        ! **** ASSEMBLING *** WRITING RANK  puts in tot matrix buffer received by idrank
                                        do ji =totistart,totiend
                                                i_contribution   = jpk*jpj_rec_a(idrank+1)*(ji-1-totistart+ relistart)
                                                do jj =totjstart,totjend
                                                        j_contribution = jpk*(jj-1-totjstart+ reljstart)
                                                        do jk =1, jpk
                                                                ind1 = jk + j_contribution + i_contribution
                                                                tottrn(jk,jj,ji)= bufftrn_TOT(ind1+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                END DO
                                unpacking_trcwrida_rank_fin_time = MPI_WTIME()

                                BeforeName = 'DA__FREQ_1/RSTbefore.'//datestring//'.'//var_to_store//'.nc'
                                BeforeNameShort = 'DA__FREQ_1/RSTbefore.'//datestring(1:11)//datestring(13:14)//datestring(16:17)//'.'//var_to_store//'.nc'
                                do ji=1,jpiglo
                                        do jj=1,jpjglo
                                                do jk=1,jpk
                                                        tottrnDA(ji,jj,jk) = REAL(tottrn(jk,jj,ji),4)
                                                end do
                                        end do
                                end do
                                writing_trcwrida_rank_init_time = MPI_WTIME()
                                CALL write_BeforeAss(BeforeName, var_to_store)
                                writing_trcwrida_rank_fin_time = MPI_Wtime()
                                !process 0 creates link to the restart file
                                !since parallel-netcdf seems to do not 
                                !read filenames with colons
                                SysErr = system("ln -sf $PWD/"//BeforeName//" "//BeforeNameShort)
                                if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
                        END IF
                        !writing_trcwrida_rank_fin_time = MPI_Wtime()

                        writing_trcwrida_rank_delta_time = writing_trcwrida_rank_fin_time - writing_trcwrida_rank_init_time
                        unpacking_trcwrida_rank_delta_time = unpacking_trcwrida_rank_fin_time - unpacking_trcwrida_rank_init_time
                        writing_trcwrida_rank_sum_time = writing_trcwrida_rank_delta_time + writing_trcwrida_rank_sum_time
                        unpacking_trcwrida_rank_sum_time = unpacking_trcwrida_rank_delta_time + unpacking_trcwrida_rank_sum_time
                        if (myrank==0) then

                                write(*,*)'writingtottime_trcwrida', writing_trcwrida_rank_sum_time,' ',jv
                                write(*,*)'unpackingtottime_trcwrida',unpacking_trcwrida_rank_sum_time,' ',jv
                        end if

                END IF
        END DO DA_LOOP
        CALL CHL_subroutine(datestring)
       
        trcwriparttime = MPI_WTIME() - trcwriparttime
        trcwritottime = trcwritottime + trcwriparttime

END SUBROUTINE trcwriDA
!------------------------------------------------------------------------------------------------
!SOLO PROC CHE FA CHL
SUBROUTINE CHL_subroutine(datestring)

        USE netcdf
        USE myalloc
        USE IO_mem
        USE calendar
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module
        USE DA_mem

        USE MPI_GATHER_INFO

        USE MATRIX_VARS
        USE NODES_MODULE
        USE DTYPE_PROCS_STRING_MODULE
        USE DA_VARS_module
        
        IMPLICIT NONE
        CHARACTER(LEN=17), INTENT(IN) :: datestring

        double precision ::  Miss_val =1.e20
        INTEGER jk,jj,ji,jn,i
        INTEGER s, nc, counter
        integer timid, depid, yid, xid, idvar
        double precision julian

        CHARACTER(LEN=45) BeforeName
        CHARACTER(LEN=43) BeforeNameShort

        CHARACTER(LEN=3) varname

        INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd,status(MPI_STATUS_SIZE)
        INTEGER irange, jrange
        INTEGER totistart, totiend, relistart, reliend
        INTEGER totjstart, totjend, reljstart, reljend
        INTEGER ind1, i_contribution, j_contribution
        INTEGER SysErr, system
        INTEGER ::jv,n_dumping_cycles,writing_rank,counter_var_DA,ivar,jn_da,ind_col
        CHARACTER(LEN=20)  var_to_store
        DOUBLE PRECISION ::chl_start,gatherv_chl_init_time,gatherv_chl_fin_time,chl_writing,chl_writing_end,chl_sum_writing,chl_end,chl_delta


        julian=datestring2sec(datestring)
        buf     = Miss_val
        bufftrn = Miss_val
        if (myrank==0) tottrn = Miss_val 

         
        chl_start = MPI_WTIME()

        if(MYRANK ==0) CHLtot = 0.0

        DA_CHL_LOOP: DO jv = 1, PX_DA
                        
                jn_da = DA_table(jv)

                !PACKING
                do ji =1 , jpi
                        i_contribution= jpk*jpj * (ji - 1 )
                        do jj =1 , jpj
                                j_contribution=jpk*(jj-1)
                                do jk =1 , jpk
                                        ind1 = jk + j_contribution + i_contribution
                                        if (tmask(jk,jj,ji).eq.1) then
                                                bufftrn(ind1)= trn(jk,jj,ji, jn_da)
                                        endif
                                enddo
                        enddo
                enddo
                gatherv_chl_init_time = MPI_Wtime()
                CALL MPI_GATHERV(bufftrn, sendcount, MPI_DOUBLE_PRECISION, bufftrn_TOT, jprcv_count, jpdispl_count, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
                gatherv_chl_fin_time = gatherv_chl_fin_time + MPI_Wtime()
                

                if(MYRANK == 0) then
                        write(*,*)'chl_gatherv_time is ',gatherv_chl_fin_time
                        var_to_store = PX_matrix(jv)%var_name

                        DO idrank = 0,mpi_glcomm_size-1

                                ! ******* WRITING RANK sets indexes of
                                ! tot matrix where to place buffers of
                                ! idrank
                                irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                totiend   = totistart + irange - 1
                                totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                totjend   = totjstart + jrange - 1
                                relistart = 1 + iPd_a(idrank+1) - 1
                                reliend   = relistart + irange - 1
                                reljstart = 1 + jPd_a(idrank+1) - 1
                                reljend   = reljstart + jrange - 1
                                ! **** ASSEMBLING *** WRITING RANK  puts
                                ! in tot matrix buffer received by
                                ! idrank
                                do ji =totistart,totiend
                                        i_contribution   = jpk*jpj_rec_a(idrank+1)*(ji-1-totistart+ relistart)
                                        do jj =totjstart,totjend
                                                j_contribution = jpk*(jj-1-totjstart+ reljstart)
                                                do jk =1, jpk
                                                        ind1 = jk + j_contribution + i_contribution
                                                        tottrn(jk,jj,ji)= bufftrn_TOT(ind1+jpdispl_count(idrank+1))
                                                enddo
                                        enddo
                                enddo
                        END DO

                        BeforeName = 'DA__FREQ_1/RSTbefore.'//datestring//'.'//var_to_store//'.nc'
                        BeforeNameShort = 'DA__FREQ_1/RSTbefore.'//datestring(1:11)//datestring(13:14)//datestring(16:17)//'.'//var_to_store//'.nc'
                        do ji=1,jpiglo
                                do jj=1,jpjglo
                                        do jk=1,jpk
                                                tottrnDA(ji,jj,jk) = REAL(tottrn(jk,jj,ji),4)
                                        end do
                                end do
                        end do
                        chl_writing = MPI_WTIME()
                        CALL write_BeforeAss(BeforeName, var_to_store)
                        chl_writing_end = MPI_WTIME()
                        chl_sum_writing = chl_sum_writing + chl_writing_end

                        write(*,*) ' chl_writing_time_is', chl_sum_writing

                        !process 0 creates link to the restart file
                        !since parallel-netcdf seems to do not 
                        !read filenames with colons
                        SysErr = system("ln -sf $PWD/"//BeforeName//" "//BeforeNameShort)
                        if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)

          
                        do jk=1,jpk
                                do jj=1,jpjglo
                                        do ji=1,jpiglo
                                                CHLtot(ji,jj,jk) = CHLtot(ji,jj,jk) + tottrnDA(ji,jj,jk)
                                        end do
                                end do
                        end do
                endif

        END DO DA_CHL_LOOP

        IF(myrank==0)then

                s = nf90_create(CHLSUP_FOR_DA, NF90_CLOBBER, nc)
                s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
                s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
                s= nf90_def_dim(nc,'depth'         , 1   ,depid)
                s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)
                s = nf90_def_var(nc, 'lchlm', nf90_float, (/xid,yid,depid,timid/),  idVAR)

                s = nf90_put_att(nc,idVAR, 'missing_value',4.e+20)

                s =nf90_enddef(nc)

                CHL_SUP = CHLtot(:,:,1)
                s = nf90_put_var(nc, idVAR  ,CHL_SUP ); call handle_err1(s,counter,CHLSUP_FOR_DA)
                s=nf90_close(nc)
        
        END IF
        chl_end = MPI_WTIME()
        chl_delta = chl_end - chl_start

        if(myrank==0)then
                write(*,*)'chl_subroutines is ', chl_delta
        end if

END SUBROUTINE CHL_subroutine

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************
!      writes tottrnDA as float on NetCDF classic file, because at the moment
!      3d_var uses parallel netcdf that does not work with NetCDF4

SUBROUTINE write_BeforeAss(fileNetCDF, VAR)

       USE netcdf
       USE myalloc
       USE DA_mem

       IMPLICIT NONE
       CHARACTER*(*), intent(in) :: fileNetCDF
       CHARACTER(LEN=3),intent(in):: VAR

       ! local

       integer s, nc, counter
       integer depid, yid, xid,idN

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'   , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'   , jpjglo,  yid)
        s= nf90_def_dim(nc,'z'   , jpk   ,depid)


        s = nf90_def_var(nc,VAR, nf90_float, (/xid,yid,depid /), idN)
        s = nf90_put_att(nc,idN   , 'missing_value',1.e+20)
        s =nf90_enddef(nc)
        s = nf90_put_var(nc, idN,  tottrnDA); call handle_err1(s,counter,fileNetCDF)
        s =nf90_close(nc)

END SUBROUTINE write_BeforeAss


