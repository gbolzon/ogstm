      SUBROUTINE trcdia(datemean, datefrom, dateend,FREQ_GROUP)

      USE myalloc
      use mpi
      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER :: IERROR
      DOUBLE PRECISION :: start_trcdit_time, end_trcdit_time,trcdit_time,trcdit_time_max, trcdit_time_sum
      DOUBLE PRECISION :: start_diadump_time, end_diadump_time, diadump_time, diadump_time_max, diadump_time_sum
       trcdiaparttime = MPI_WTIME() ! cronometer-start

!      writes ave files for tracer concentration

      start_trcdit_time = MPI_Wtime()
      CALL   trcdit(datemean, datefrom, dateend,FREQ_GROUP)
      end_trcdit_time = MPI_Wtime()
      trcdit_time = end_trcdit_time - start_trcdit_time 
      CALL MPI_Reduce(trcdit_time,trcdit_time_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_Reduce(trcdit_time,trcdit_time_sum,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
      if (myrank==0) then
                write(*,*)'TRCDIT_MAX_TIME is',trcdit_time_max
                write(*,*)'TRCDIT_MEAN_TIME is',trcdit_time_sum/mpi_glcomm_size
      end if

      start_diadump_time = MPI_Wtime()
      CALL  diadump(datemean, datefrom, dateend,FREQ_GROUP)
      end_diadump_time = MPI_Wtime()
      diadump_time = end_diadump_time - start_diadump_time    
      CALL MPI_Reduce(diadump_time,diadump_time_max,3,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_Reduce(diadump_time,diadump_time_sum,4,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
      if (myrank==0) then
                write(*,*)'DIADUMP_MAX_TIME is',diadump_time_max
                write(*,*)'DIADUMP_MEAN_TIME is',diadump_time_sum/mpi_glcomm_size
      end if

      CALL fluxdump(datemean, datefrom, dateend,FREQ_GROUP)


       trcdiaparttime =   MPI_WTIME() - trcdiaparttime  ! cronometer-stop
       trcdiatottime  = trcdiatottime + trcdiaparttime

      END SUBROUTINE trcdia
