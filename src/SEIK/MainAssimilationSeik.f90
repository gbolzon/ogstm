      SUBROUTINE MainAssimilationSeik(datestr)
      use myalloc
      use filenames
      use DA_mem
      use mpi
      !use DA_Params
      !use ogstm_mpi_module
      !use mpi_str, only: Var3DCommunicator
      !use calendar

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestr

      character(LEN=1024) SATFILE, VARFILE
      character(LEN=2) MONTH
      character (LEN=8) DAY
      character MISFIT_OPT
      logical ISLOG
      integer :: ierr !, ii
      !integer :: SysErr, system
        
        INTEGER jk,jj,ji,jn
        double precision :: fillValue, fillvalue999

      DAparttime  = MPI_WTIME()
      MONTH=datestr(5:6)
      DAY  =datestr(1:8)

      ! DA_Params init
      !DA_Date = datestr
      !ShortDate = datestr(1:11)//datestr(13:14)//datestr(16:17)
      !jpk_200 = AssimilationLevels
      !NBioVar = 17
      !DA_JulianDate = datestring2sec(datestr)

      !MISFIT_OPT ='2'
      !ISLOG = .false.
      
      MISFIT_OPT ='1'
      ISLOG = .true.

      ! SATDIR li vuole pronti all'uso, già tagliati e interpolati
      SATFILE   = 'SATELLITE/' // DAY // trim(satfile_suffix)
      VARFILE   = 'DA_static_data/MISFIT/VAR2D/var2D.' // MONTH // '.nc'
      !EOF_FILE  = 'DA_static_data/3D_VAR/EOF/eof.'  // MONTH // '.nc'
      !GRID_FILE = 'DA_static_data/3D_VAR/GRID/BFM_grid.nc'
      !ANIS_FILE = 'DA_static_data/3D_VAR/gradsal.nc'

      MISFIT_FILE='DA__FREQ_1/'// DAY // '.chl_mis.nc'
      !ARGO_FILE='DA__FREQ_1/'// DAY // '.P_l_arg_mis.dat'
      !CORR_FILE = 'DA__FREQ_1/'// DAY // '_corr.nc'
      !EIV_FILE  = 'DA__FREQ_1/'// DAY // '_eiv.nc'
      !OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'


      CHLSUP_FOR_DA = 'DA__FREQ_1/chl.' // datestr // '.nc'
      if (EnsembleRank==0) then

!write(*,*) EnsembleRank, myrank, "analisi0"

        CALL trcwriDA(DATEstr)  ! Dumps Before Assimilation real*4

!call mpi_barrier(LocalComm, ierr)
!write(*,*) EnsembleRank, myrank, "analisi1"	        
!call mpi_barrier(LocalComm, ierr)
!call mpi_abort(mpi_comm_world, -1, ierr)

      !if (myrank .lt. DA_Nprocs ) then

          !allocate(DA_VarList(NBioVar))
          !do ii=1,NBioVar
          !  DA_VarList(ii) = varlistDA(ii)
          !enddo

          if(myrank .eq. 0) then ! .and. sat .eq. 1) then
            write(*,*) 'satfile=', trim(SATFILE)
            write(*,*) 'varfile=', trim(VARFILE)
            write(*,*) 'misfit=', trim(MISFIT_FILE)
            call CREATEMISFIT(SATFILE,VARFILE,MISFIT_OPT, ISLOG, MISFIT_FILE) ! produces MISFIT.nc
!write(*,*) EnsembleRank, myrank, "analisi2"
            !write(*,*) 'eof = ',   trim(EOF_FILE)
            !write(*,*) 'grid = ',  trim(GRID_FILE)
          endif

          ! if(myrank .eq. 0) then
          !   SysErr = system("../float_preproc/Float_misfit_gen.sh -d ../float_preproc -t "//DAY)
          !   if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
          ! endif
  
          !call MPI_Barrier(Var3DCommunicator, ierr)


          !call BIOVAR

          ! deallocation of DA_VarList array
          !call clean_da_params
      endif

!call mpi_barrier(mpi_comm_world, ierr)
!write(*,*) EnsembleRank, myrank, "analisi3"
!call mpi_barrier(LocalComm, ierr)

!----------New Seik Part--------------------------------------------------------
      
if (SeikDim>0) then
        ComputedObsSeik=0.0d0
        MisfitSeik=0.0d0
        do jn=1, jptra
            if (isaCHLVAR(ctrcnm(jn))) then
                ComputedObsSeik=ComputedObsSeik+trnEnsemble(1,:,:,jn)
                MisfitSeik=MisfitSeik+trn(1,:,:,jn)
            end if
        end do
        !ComputedObsSeik=ComputedObsSeik*bfmmask(1,:,:)
        
        !ObsBaseMember=1.0d0/log(1.3d0)**2
        !call CutCellsSurface(ObsBaseMember)
        !ObsBaseMember=ObsBaseMember*SeikMask(1,:,:)
        !ObsErrorDiag1=reshape(ObsBaseMember,(/ObsSpaceDim/))

        call readnc_slice_double_2d(trim(SATFILE),trim(satvarname), ObsDataSeik)

!call mpi_barrier(mpi_comm_world, ierr)
!write(*,*) EnsembleRank, myrank, "analisi4"
!call mpi_barrier(LocalComm, ierr)
!call mpi_abort(mpi_comm_world, -1, ierr)

        fillvalue999=-999.0d0
        fillValue = 1.0d20

        do ji=1,jpi
            do jj=1,jpj
                if ( isnanSeik(ObsDataSeik(jj,ji)).or.(ObsDataSeik(jj,ji).eq.fillValue).or. & 
                    (ObsDataSeik(jj,ji).eq.fillvalue999).or.(SeikMask(1,jj,ji).eq.0)) then !SeikMask instead of tmask instead of bfmmask
                    ComputedObsSeik(jj,ji)=0.0d0
                    ObsDataSeik(jj,ji)=0.0d0
                    MisfitSeik(jj,ji)=0.0d0
                else
                    if (ComputedObsSeik(jj,ji)>0.01d0) then
                        ComputedObsSeik(jj,ji)=log(ComputedObsSeik(jj,ji))
                    else
                        ComputedObsSeik(jj,ji)=log(0.01d0)
                    end if
                    
                    if (ObsDataSeik(jj,ji)>0.01d0) then
                        ObsDataSeik(jj,ji)=log(ObsDataSeik(jj,ji))
                    else
                        ObsDataSeik(jj,ji)=log(0.01d0)
                    end if
                    
                    if (MisfitSeik(jj,ji)>0.01d0) then
                        MisfitSeik(jj,ji)=log(MisfitSeik(jj,ji))
                    else
                        MisfitSeik(jj,ji)=log(0.01d0)
                    end if
                    MisfitSeik(jj,ji)=ObsDataSeik(jj,ji)-MisfitSeik(jj,ji)
                end if
            end do
        end do

!call mpi_barrier(mpi_comm_world, ierr)
!write(*,*) EnsembleRank, myrank, "analisi5"
!call mpi_barrier(LocalComm, ierr)
        
        if (UseModSeik) then
            call ModSeikAnalysis
        else
            call SeikAnalysis
        end if
        
!call mpi_barrier(mpi_comm_world, ierr)
!write(*,*) EnsembleRank, myrank, "analisi6"
!call mpi_barrier(LocalComm, ierr)
        
        if (EnsembleRank==0) then
            call trcwriSeik(DATEstr, -1, 'DA_SEIK/', trn)
            if (myrank==0) call WriteCov1Seik('DA_SEIK/COV1.'//DateStr//'.csv')
            if (SeikDim==1) call trcwriSeik(DATEstr, -1, 'DA_SEIK/RATIO/', trnEnsemble) ! only if we dont' have a 3rd ensemble member for parallel writing
        else if (EnsembleRank==1) then
            trb=trn-trb
            call trcwriSeik(DATEstr, -1, 'DA_SEIK/DIFFERENCE/', trb)
        else if (EnsembleRank==2) then
            call trcwriSeik(DATEstr, -1, 'DA_SEIK/RATIO/', trnEnsemble)
        end if
        
        trb=trn

end if

      !call mppsync()
      !CALL trcrstDA(datestr)

      DAparttime  = MPI_WTIME() - DAparttime
      DAtottime   = DAtottime + DAparttime
      
      ! Sviluppo : se non serve il RST After Assimilation, puo' fare uno scatter della nuova variabile di stato      
      CONTAINS

        ! ***********************************************************************
        LOGICAL FUNCTION isnanSeik(A)
        implicit none
        double precision, INTENT(IN) :: A
        if ( A.eq.A ) then
          isnanSeik = .FALSE.
          else
          isnanSeik = .TRUE.
        end if
        END FUNCTION isnanSeik
      END SUBROUTINE MainAssimilationSeik
 
