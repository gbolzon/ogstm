MODULE module_step

      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE BC_mem
      USE IO_mem, only: ave_counter_1, ave_counter_2
      USE mpi
      USE mod_atmbc
      USE mod_cbc
      ! USE mod_gibbc
      ! USE mod_tinbc

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      USE ogstm_mpi_module
use StringSEIK

 implicit NONE

 contains
 
 SUBROUTINE step
!---------------------------------------------------------------------
!
!                       ROUTINE STEP
!                     ****************
!
!  PURPOSE :
!  ---------
!    Time loop of ogstm
!
!   METHOD :
!   -------
!
!      READ the dynamics fiels and update for time step
!
!      Call passive tracer model trcstp
!                     and diagnostics trcdia
!                     and specific restart trcwri
!
!   INPUT :
!   -----
!
!
!   OUTPUT :            : no
!   ------
!
!   EXTERNAL :

!      trcstp, trcdia  passive tracers interface

       IMPLICIT NONE


! local declarations
! ==================
      INTEGER TAU, indic

      character(LEN=17)  datestring, datemean, datefrom_1, datefrom_2, NextDateString
      character(LEN=17)  date_aveforDA
      double precision sec
      LOGICAL B, isFIRST
      INTEGER :: jk,jj,ji,jn
      integer ierr
!++++++++++++++++++++++++++++++c
!         Time  loop           c
!++++++++++++++++++++++++++++++c

#ifdef ExecDA
!if (PCANeeded.ne."read") then
#endif

      isFIRST=.true.

      indic=1  ! 1 for indic to initialize output files (first call)

      TauAVEfrom_1 = TimeStepStart 
       if (IsStartBackup_1) TauAVEfrom_1 = datestringToTAU(BKPdatefrom_1)
      TauAVEfrom_2 = TimeStepStart 
       if (IsStartBackup_2) TauAVEfrom_2 = datestringToTAU(BKPdatefrom_2)
       
!call mpi_barrier(mpi_comm_world,ierr)
!where (trn<1.0d-10) trn=1.0d-10
!trb=trn

      DO TAU = TimeStepStart, TimeStep__End

         stpparttime = MPI_WTIME()  ! stop cronomether
         call tau2datestring(TAU, DATEstring)
         call tau2datestring(TAU+1,NextDateString)
         sec=datestring2sec(DATEstring)

         NOW_datestring = DATEstring ! update time manager module
         NOW_sec        = sec
         COMMON_DATESTRING = DATEstring

         call yearly(DATEstring) ! Performs yearly updates
         call daily(DATEstring)  ! Performs daily  updates


         if(lwp) write(numout,'(A,I8,A,A)') "step ------------ Starting timestep = ",TAU,' time ',DATEstring
         if(lwp) write(*,'(A,I8,A,A)')      "step ------------ Starting timestep = ",TAU,' time ',DATEstring
         !write(*,*) is_night(DATEstring)
         
#ifdef ExecDA
        if ((PCANeeded.eq."write").and.(EnsembleRank==0)) call PCACreateMatrices(DATEstring)

        if (IsaRestart(DATEstring).and.(SeikDim.gt.0))  then
            call WriteBaseSeik(datestring)
            if (EnsembleRank==0) call trcwriSeik(dateString, -1, "REDUCED_BASE/DIMENSION_"//int2str(SeikDim,3) //"/", trnVariance) !questa e' la varianza del giorno prima
        endif
#endif

        if (IsaRestart(DATEstring).and.(EnsembleRank==0))  then
            !CALL trcwri(DATEstring) ! writes the restart files
            call trcwriSeik(DATEstring, -1, 'RESTARTS/', trn)

         
            if (AVE_FREQ1%N .gt.0) then              !  void 1.aveTimes -> no backup
            if (.not.IsAnAveDump(DATEstring,1)) then ! backup conditions group 1
               call tau2datestring(TauAVEfrom_1, datefrom_1)
               CALL trcdia(datestring, datefrom_1, datestring,1)
            endif
            endif

            if (AVE_FREQ2%N .gt.0) then
            if (.not.IsAnAveDump(DATEstring,2)) then ! backup conditions group 2
               call tau2datestring(TauAVEfrom_2, datefrom_2)
               if (save_bkp_group2) CALL trcdia(datestring, datefrom_2, datestring,2)
            endif
            endif

         if (lwp) then
             B = writeTemporization("trcdia____", trcdiatottime)
             B = writeTemporization("trcwri____", trcwritottime)
         endif
         endif



! For offline simulation READ DATA or precalculalted dynamics fields
! ------------------------------------------------------------------

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "written"

      CALL forcings_PHYS(DATEstring)

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "forcings_PHYS"

      CALL forcings_KEXT(datestring)

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "forcings_KEXT"

      !CALL bc_gib       (DATEstring)     ! CALL dtatrc(istp,0)! Gibraltar strait BC
      !CALL bc_tin       (DATEstring)     ! CALL dtatrc(istp,1)

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      !bc_tin_partTime = MPI_WTIME()
      call boundaries%update(datestring)

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "boundaries"

      !bc_tin_partTime = MPI_WTIME()    - bc_tin_partTime
      !bc_tin_TotTime  = bc_tin_TotTime + bc_tin_partTime

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      CALL bc_atm       (DATEstring)     ! CALL dtatrc(istp,2)

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "bc_atm"

      CALL bc_co2       (DATEstring)

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "bc_co2"

      CALL eos          ()               ! Water density

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "eos"

      if (IsAnAveDump(DATEstring,1).and.(EnsembleRank==0)) then
         call MIDDLEDATE(TauAVEfrom_1, TAU, datemean)

         call tau2datestring(TauAVEfrom_1, datefrom_1)
         if (IsStartBackup_1) datefrom_1 = BKPdatefrom_1 ! overwrite
         CALL trcdia(datemean, datefrom_1, datestring,1)
         TauAVEfrom_1    = TAU
         ave_counter_1   = 0   !  reset the counter
         IsStartBackup_1 = .false.
        if (lwp)  B = writeTemporization("trcdia____", trcdiatottime)
      endif


      if (IsAnAveDump(DATEstring,2).and.(EnsembleRank==0)) then
         call MIDDLEDATE(TauAVEfrom_2, TAU, datemean)

         call tau2datestring(TauAVEfrom_2, datefrom_2)
         if (IsStartBackup_2) datefrom_2 = BKPdatefrom_2 ! overwrite
         CALL trcdia(datemean, datefrom_2, datestring,2)
         TauAVEfrom_2    = TAU
         ave_counter_2   = 0   !  reset the counter
         IsStartBackup_2 = .false.
         if (lwp) B = writeTemporization("trcdia____", trcdiatottime)
      endif



#ifdef ExecDA_I_dont_want_this
! 3D-VAR Assimilation
    if (SeikDim.eq.0) then 
      if (IsaDataAssimilation(DATEstring)) then
        call tau2datestring(TauAVEfrom_1, datefrom_1)
        CALL mainAssimilation(DATEstring, datefrom_1)
         if (lwp) B = writeTemporization("DATA_ASSIMILATION____", DAparttime)
      endif
    end if
#endif

#ifdef ExecDA

! Seik Assimilation

        !if ((SeikDim.gt.0).and.(IsaDataAssimilation(DATEstring))) then            
        if (IsaDataAssimilation(DATEstring)) then
            CALL MainAssimilationSeik(DATEstring)            
            if (lwp) B = writeTemporization("DATA_ASSIMILATION____", DAparttime)
        endif

! Ensemble creation

        if ((SeikDim.gt.0).and.(datestring(10:17).eq."00:00:00")) then
            call SeikCreateEnsemble()

            if (.false.) then !.true. if u want to save the ensemble before the evolution
                call trcwriSeik(datestring, EnsembleRank, 'ENSEMBLE/', trn)
            endif

        endif
#endif

!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "ensemble"

if (.false.) then
trnEnsemble=trn
!if(lwp) write(*,*) ctrcnm
do jn=1, jptra
do ji=1, jpi
do jj=1,jpj
do jk=1,jpk
    if ((trn(jk,jj,ji,jn)<1.0d-8).and.(bfmmask(jk,jj,ji)==1)) then
        write(*,*) "Small:en", EnsembleRank, "my" , MyRank, "n ", ctrcnm(jn), "(",jk,jj,ji,jn,")=",trn(:,jj,ji,jn) 
        !write(*,*) trn(jk,jj,ji,:)
        exit
    end if
end do
end do
end do
end do
end if

!where ((trn<1.0d-12).or.(.not.(trn==trn))) trn=1.0d-12
!trb=trn

! Call Passive tracer model between synchronization for small parallelisation
        CALL trcstp    ! se commento questo non fa calcoli
        
if (.false.) then
do jn=1, jptra
do ji=1,jpi
do jj=1,jpj
do jk=1,jpk
    if ((trn(jk,jj,ji,jn)<trnEnsemble(jk,jj,ji,jn)*0.5).and.(bfmmask(jk,jj,ji)==1)) then
        write(*,*) "Halved:en", EnsembleRank, "my" , MyRank, "n ", ctrcnm(jn), "(",jk,jj,ji,jn,")=",trn(:,jj,ji,jn)
        !write(*,*) trn(jk,jj,ji,:)
        exit
    end if
end do
end do
end do
end do
end if

#ifdef ExecDA
        if ((SeikDim.gt.0).and.(NextDateString(10:17).eq."00:00:00")) then

            if (.false.) then !.true. if you want to save the ensamble after the evolution
                call trcwriSeik(datestring, EnsembleRank, 'PENSEMBLE/', trn)
            endif

            call SeikForecast(datestring)

        endif
#endif


        if (EnsembleRank==0) call trcave
        ave_counter_1 = ave_counter_1 +1  ! incrementing our counters
        ave_counter_2 = ave_counter_2 +1

call mpi_barrier(mpi_comm_world, ierr)
       stpparttime = MPI_WTIME() - stpparttime
       stptottime  = stptottime  + stpparttime
        if (lwp) write(*,*) "Step in sec ", stpparttime


! OGSTM TEMPORIZATION
       IF (TAU.GT.TimeStepStart) THEN
        IF( mod( TAU, nwritetrc ).EQ.0) THEN
           if (lwp) then
               write(*,*) "************* OGSTM TEMPORIZATION *********"
               write(*,*) "              Iteration",TAU
               write(*,*) "routine******time_tot*********time_ave*****"
           endif
           B = writeTemporization("forPhys___", forcing_phys_TotTime)
           B = writeTemporization("forKext___", forcing_kext_TotTime)
           B = writeTemporization("bcCO2_____", bc_co2_TotTime)
           B = writeTemporization("bcTIN_____", bc_tin_TotTime)
           B = writeTemporization("bcATM_____", bc_atm_TotTime)
           B = writeTemporization("bcGIB_____", bc_gib_TotTime)
           B = writeTemporization("density___", density_TotTime)
           B = writeTemporization("averaging_", ave_TotTime   )
           B = writeTemporization("trcopt____", trcopttottime)
           B = writeTemporization("trcbio____", BIOtottime)
           B = writeTemporization("trcadv____", trcadvtottime)
           B = writeTemporization("trcdmp____", trcdmptottime)
           B = writeTemporization("trcbil____", trcbilaphdftottime)
           B = writeTemporization("trcsbc____", trcsbctottime)
           B = writeTemporization("trcsms____", trcsmstottime)
           B = writeTemporization("trczdf____", trczdftottime)
           B = writeTemporization("snutel____", snuteltottime)
           B = writeTemporization("check_____", checkVtottime)
           B = writeTemporization("trcnxt____", trcnxttottime)
           B = writeTemporization("trcstp____", trcstptottime)


           B = writeTemporization("flxdump___",flx_TotTime  )
           B = writeTemporization("stp_______", stptottime  )

           call reset_Timers()
       ENDIF
      ENDIF


!+++++++++++++++++++++++++++++c
!      End of time loop       c
!+++++++++++++++++++++++++++++c

      END DO  
      
! PCA if needed
#ifdef ExecDA
!end if
        if ((PCANeeded.eq."read").and.(EnsembleRank==0)) call PCAReadMatrices
        if (((PCANeeded.eq."read").or.(PCANeeded.eq."write")).and.(EnsembleRank==0)) call PCASeik
#endif

      CONTAINS

      LOGICAL FUNCTION writeTemporization(string, elapsedtime)
      IMPLICIT NONE
      CHARACTER(LEN=*) string
      double precision elapsedtime

      if (isFIRST) then
         write(*,250) string,elapsedtime,elapsedtime/(TAU-TimeStepStart +1)," myrank->", myrank
         isFirst=.false.
      else
         write(*,250) string,elapsedtime,elapsedtime/nwritetrc," myrank->", myrank
      endif
      writeTemporization = .true.
250   FORMAT (A , ES11.4 ,ES20.7 ,A20 , I3 )
      END FUNCTION writeTemporization

      END SUBROUTINE step


      SUBROUTINE trcstp
!---------------------------------------------------------------------
!
!                       ROUTINE trcstp
!                     *****************
!
!  PURPOSE :
!  ---------
!	time loop of ogstm for passive tracer
!
!   METHOD :
!   -------
!      compute the well/spring evolution
!
!      compute the time evolution of tracers concentration
!         with advection
!         with horizontal diffusion
!         with surface boundary condition
!         with IMPLICIT vertical diffusion

       IMPLICIT NONE
      integer jn,jk,ji,jj
integer ierr

      trcstpparttime = MPI_WTIME() ! cronometer-start

      IF (ladv) CALL trcadv ! tracers: advection

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "adv"

#    if defined key_trc_dmp
      CALL trcdmp ! tracers: damping for passive tracerstrcstp

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "dmp"


! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call boundaries%apply(e3t, trb, tra)

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "bound"


! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

#    endif

! tracers: horizontal diffusion IF namelist flags are activated
! -----------------------------

      IF (lhdf) CALL trchdf

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "hdf"


! tracers: sink and source (must be  parallelized on vertical slab)
      IF (lsbc) CALL trcsbc ! surface cell processes, default lsbc = False

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "sbc"


      IF (lbfm) CALL trcsms

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "sms"


      IF (lzdf) CALL trczdf ! tracers: vertical diffusion

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "zdf"


      IF (lsnu) CALL snutel

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "snute"


      IF (lhtp) CALL hard_tissue_pump

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "hard"


      ! CALL checkValues

      CALL trcnxt ! tracers: fields at next time step

!call prova
!call mpi_barrier(mpi_comm_world,ierr)
!if (lwp) write(*,*) "nxt"

      
      trcstpparttime = MPI_WTIME() - trcstpparttime ! cronometer-stop
      trcstptottime = trcstptottime + trcstpparttime

      END SUBROUTINE trcstp

subroutine prova
implicit none
integer ji,jj,jk,jn

if(.false.) then
do jn=1, jptra
do ji=1,jpi
do jj=1,jpj
do jk=1,jpk
    if ((trn(jk,jj,ji,jn)<trnEnsemble(jk,jj,ji,jn)*0.5).and.(bfmmask(jk,jj,ji)==1)) then
        write(*,*) "Halved:en", EnsembleRank, "my" , MyRank, "n ", ctrcnm(jn), "(",jk,jj,ji,jn,")=",trn(:,jj,ji,jn)
        !write(*,*) trn(jk,jj,ji,:)
        exit
    end if
end do
end do
end do
end do
end if

end subroutine



end module
