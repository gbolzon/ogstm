      SUBROUTINE MainAssmimilationSeik(datestr)
      use myalloc
      use filenames
      use DA_mem
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
      !integer :: ierr, ii
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
        CALL trcwriDA(DATEstr)  ! Dumps Before Assimilation real*4

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
      
!----------New Seik Part--------------------------------------------------------
      
        ComputedObsSeik=0.0d0
        MisfitSeik=0.0d0
        do jn=1, jptra
            if (isaCHLVAR(ctrcnm(jn))) then
                ComputedObsSeik=ComputedObsSeik+trnEnsemble(1,:,:,jn)
                MisfitSeik=MisfitSeik+trn(1,:,:,jn)
        end do
        !ComputedObsSeik=ComputedObsSeik*bfmmask(1,:,:)
        
        ObsErrorDiag1=1.0d0/log(1.1d0)**2
        !ObsErrorDiag1=1.0d0
        
        call readnc_slice_double_2d(trim(SATFILE),trim(satvarname), ObsDataSeik)
        fillvalue999=-999.0d0
        fillValue = 1.0d20

        do ji=1,jpi
            do jj=1,jpj
                if ( isnan2(ObsDataSeik(jj,ji)).or.(ObsDataSeik(jj,ji).eq.fillValue).or. & 
                    (ObsDataSeik(jj,ji).eq.fillvalue999).or.(bfmmask(1,jj,ji).eq.0)) then
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
                        MisfitSeik(jj,ji)=ObsDataSeik(jj,ji)-MisfitSeik(jj,ji)
                    else
                        MisfitSeik(jj,ji)=log(0.01d0)
                    end if
                end if
            end do
        end do

        call SeikAnalysis
        
        if (EnsembleRank==0) then
            call trcwriSeik(DATEstring, -1, 'DA_SEIK/', trn)
            if (myrank==0) call WriteCov1Seik('DA_SEIK/COV1.'//DateString//'.csv')
            if (SeikDim==1) trcwriSeik(DATEstring, -1, 'DA_SEIK/RATIO/', Basemember) ! only if we dont' have a 3rd ensemble member for parallel writing
        else if (EnsembleRank==1) then
            trb=trn-trb
            call trcwriSeik(DATEstring, -1, 'DA_SEIK/DIFFERENCE/', trn)
        else if (EnsembleRank==2) then
            call trcwriSeik(DATEstring, -1, 'DA_SEIK/RATIO/', Basemember)
        end if
        
        trb=trn

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
          isnan2 = .FALSE.
          else
          isnan2 = .TRUE.
        end if
        END FUNCTION isnanSeik
      END SUBROUTINE MainAssmimilationSeik
 
