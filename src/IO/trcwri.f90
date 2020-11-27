        SUBROUTINE trcwri(datestring)

        USE myalloc
        USE IO_mem
        USE calendar
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module

        USE MPI_GATHER_INFO

        USE MATRIX_VARS
        USE NODES_MODULE
        USE DTYPE_PROCS_STRING_MODULE



        IMPLICIT NONE
        CHARACTER(LEN=17), INTENT(IN) :: datestring

!----------------------------------------------------------------------
! local declarations
! ==================
        double precision ::  Miss_val =1.e20
        INTEGER jk,jj,ji,i
        double precision julian


        CHARACTER(LEN=37) filename

        CHARACTER(LEN=3) varname

        INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
        INTEGER irange, jrange
        INTEGER totistart, totiend, relistart, reliend
        INTEGER totjstart, totjend, reljstart, reljend
        INTEGER ind1, i_contribution, j_contribution
        CHARACTER(LEN=20)  var_to_store
        INTEGER :: COUNTER_VAR_TRCWRI, n_dumping_cycles, jv, ivar, writing_rank, ind_col 
        DOUBLE PRECISION :: trcwri_packing_init_time, start_time_trcwri, packing_init_time_trcwri, trcwri_gatherv_init_time, trcwri_gatherv_fin_time, trcwri_packing_fin_time, trcwri_gatherv_delta_time, trcwri_packing_delta_time, trcwri_gatherv_sum_time, trcwri_packing_sum_time, trcwri_writing_rank_init_time, trcwri_unpacking_rank_init_time, trcwri_unpacking_rank_fin_time, trcwri_writing_rank_fin_time, trcwri_writing_rank_delta_time, trcwri_unpacking_rank_delta_time
        DOUBLE PRECISION :: trcwri_gatherv_cycles, trcwri_writing_rank_sum_time, trcwri_unpacking_rank_sum_time, trcwri_gatherv_delta1_time

        filename = 'RST.20111231-15:30:00.N1p.nc'
        julian=datestring2sec(datestring)

        trcwriparttime = MPI_WTIME() ! cronometer-start

        call mppsync()
        start_time_trcwri =  MPI_Wtime()
        buf = Miss_val
        bufftrn = Miss_val
        !trcwri_gatherv_fin_time = 0
        trcwri_gatherv_cycles = 0
        trcwri_writing_rank_sum_time = 0
        trcwri_unpacking_rank_sum_time = 0
        if (WRITING_RANK_WR) tottrn = Miss_val


        n_dumping_cycles = matrix_state_2_row

        COUNTER_VAR_TRCWRI = 1

        RESTARTS_LOOP: DO jv = 1, n_dumping_cycles
                                
                DO ivar = 1 , nodes*num_of_wr_procs_perNODE
        
                        writing_rank = writing_procs(ivar)
                        IF (COUNTER_VAR_TRCWRI > JPTRA)then
                                EXIT
                        ELSE

                                do ji =1 , jpi
                                        i_contribution= jpk*jpj * (ji - 1 )
                                        do jj =1 , jpj
                                                j_contribution=jpk*(jj-1)
                                                do jk =1 , jpk
                                                        ind1 = jk + j_contribution + i_contribution
                                                        if (tmask(jk,jj,ji).eq.1) then
                                                                bufftrn(ind1)= trn(jk,jj,ji,COUNTER_VAR_TRCWRI)
                                                        endif
                                                enddo
                                        enddo
                                enddo
                                trcwri_packing_fin_time = MPI_Wtime()
                                trcwri_packing_delta_time = trcwri_packing_fin_time - trcwri_packing_init_time
                                
                                COUNTER_VAR_TRCWRI = COUNTER_VAR_TRCWRI + 1
                                trcwri_gatherv_init_time = MPI_Wtime()
                                CALL MPI_GATHERV(bufftrn, sendcount, MPI_DOUBLE_PRECISION, bufftrn_TOT,jprcv_count, jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,MPI_COMM_WORLD, IERR)
                                trcwri_gatherv_fin_time = MPI_Wtime()
                                trcwri_gatherv_delta_time = trcwri_gatherv_fin_time - trcwri_gatherv_init_time
                                !trcwri_gatherv_cycles = trcwri_gatherv_delta1_time + trcwri_gatherv_cycles
                                if(myrank==0)then
                                        write(*,*)'trcwrigathervDElt',trcwri_gatherv_delta_time
                                        write(*,*)'trcwri_packing',trcwri_packing_delta_time
                                end if
                        END IF

                END DO
                !trcwri_packing_fin_time = MPI_Wtime()

                !trcwri_gatherv_delta_time = trcwri_gatherv_cycles - trcwri_gatherv_init_time
                !trcwri_packing_delta_time = trcwri_packing_fin_time - trcwri_packing_init_time
                !CALL MPI_Reduce( trcwri_gatherv_delta_time, trcwri_gatherv_sum_time,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,IERROR)
                !CALL MPI_Reduce( trcwri_packing_delta_time, trcwri_packing_sum_time,2,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,IERROR)
                !if(myrank==0)then
                !        write(*,*)'trcwri_gatherv_sum_time_trcdit is ',trcwri_gatherv_sum_time
                !        write(*,*)'trcwri_packing_sum_time_trcdit is ',trcwri_packing_sum_time
                !end if


                if(WRITING_RANK_WR) then
                        !trcwri_writing_rank_init_time = MPI_Wtime()
                        trcwri_unpacking_rank_init_time = MPI_WTIME()
                        do i=1, nodes*num_of_wr_procs_perNODE
                                if(myrank == writing_procs(i))then
                                        ind_col=i
                                        exit
                                end if
                        end do
                        !ind_col = (myrank / n_ranks_per_node)+1
                        var_to_store =matrix_state_2(jv,ind_col)%var_name
                        IF (var_to_store == "novars_input")then
                                EXIT
                        ELSE
                                DO idrank = 0,mpi_glcomm_size-1
                                        ! ******* WRITING RANK sets
                                        ! indexes of tot matrix where to
                                        ! place buffers of idrank
                                        irange    = iPe_a(idrank+1) -iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) -jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) +iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange -1
                                        totjstart = jstart_a(idrank+1) +jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange -1
                                        relistart = 1 + iPd_a(idrank+1)- 1
                                        reliend   = relistart + irange -1
                                        reljstart = 1 + jPd_a(idrank+1)- 1
                                        reljend   = reljstart + jrange -1
                                        ! **** ASSEMBLING *** WRITING
                                        ! RANK  puts in tot matrix
                                        ! buffer received by idrank
                                        do ji =totistart,totiend
                                                i_contribution   =jpk*jpj_rec_a(idrank+1)*(ji-1-totistart+ relistart)
                                                do jj =totjstart,totjend
                                                        j_contribution =jpk*(jj-1-totjstart+ reljstart)
                                                        do jk =1, jpk
                                                                ind1 =jk + j_contribution + i_contribution
                                                                tottrn(jk,jj,ji)=bufftrn_TOT(ind1+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                END DO
                                trcwri_unpacking_rank_fin_time = MPI_WTIME()
                                filename = 'RESTARTS/RST.'//datestring//'.'//trim(var_to_store)//'.nc'
                                trcwri_unpacking_rank_delta_time = trcwri_unpacking_rank_fin_time - trcwri_unpacking_rank_init_time
                                trcwri_writing_rank_init_time = MPI_Wtime()
                                CALL write_restart(filename,var_to_store,julian, deflate_rst, deflate_level_rst)
                                trcwri_writing_rank_fin_time = MPI_Wtime()
                                trcwri_writing_rank_delta_time = trcwri_writing_rank_fin_time - trcwri_writing_rank_init_time
                                if (myrank==0) then

                                        write(*,*)'trcwri_writingtottime',trcwri_writing_rank_delta_time,'',jv
                                        write(*,*)'trcwri_unpackingtottime',trcwri_unpacking_rank_delta_time,' ',jv
                                end if
       
                        END IF
                        !trcwri_writing_rank_fin_time = MPI_Wtime()

                        !trcwri_writing_rank_delta_time = trcwri_writing_rank_fin_time - trcwri_writing_rank_init_time
                        !trcwri_unpacking_rank_delta_time = trcwri_unpacking_rank_fin_time - trcwri_unpacking_rank_init_time
                        !trcwri_writing_rank_sum_time = trcwri_writing_rank_delta_time + trcwri_writing_rank_sum_time
                        !trcwri_unpacking_rank_sum_time = trcwri_unpacking_rank_delta_time + trcwri_unpacking_rank_sum_time
                        !if (myrank==0) then
!
!                                write(*,*)'trcwri_writingtottime',trcwri_writing_rank_sum_time,' ',jv
!                                write(*,*)'trcwri_unpackingtottime',trcwri_unpacking_rank_sum_time, ' ',jv
!                        end if
                END IF
        END DO RESTARTS_LOOP
        END SUBROUTINE trcwri
