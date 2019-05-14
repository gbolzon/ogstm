SUBROUTINE ReadBaseSeik
!---------------------------------------------------------------------
!
!                       ROUTINE ReadBaseSeik
!                     ************************
!
!  PURPOSE :
!  ---------
!     READ files of SEIK's reduced base 
!
!----------------------------------------------------------------------


    USE calendar
    USE myalloc
    USE TIME_MANAGER
    use mpi

    IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
    INTEGER :: jn, BaseIndex, ierr
    CHARACTER(LEN=26) :: DirName
    CHARACTER(LEN=59) :: FileNameBase
    CHARACTER(LEN=53) :: FileNameCov

    DirName='REDUCED_BASE/DIMENSION_010'
    FileNameBase='REDUCED_BASE/DIMENSION_010/BASE001.20111231-15:30:00.N1p.nc'
    FileNameCov='REDUCED_BASE/DIMENSION_010/COV1.20111231-15:30:00.csv'
    
    write(DirName,'(A23,I3.3)') 'REDUCED_BASE/DIMENSION_', SeikDim
    if (EnsembleRank==NotWorkingMember) then
        if (myrank==0) then
            FileNameCov=DirName//'/COV1.'//DateStart//'.csv'
            open(unit=CovIOUnit, file=FileNameCov, form='formatted', iostat=ierr, action='read', access='sequential',status='old')
            if (ierr/=0) then
                write(*,*) 'Error initializing covariance matrix from file: ', ierr
                write(*,*) 'I will stop'
                call mpi_abort(mpi_comm_world,1,ierr)
            end if
            do jn=1,SeikDim
                read(CovIOUnit,*,iostat=ierr) CovSeik1(:,jn)
                if (ierr/=0) then
                    write(*,*) 'Error initializing covariance matrix from file, while reading line', jn, 'error:', ierr
                    write(*,*) 'I will stop'
                    call mpi_abort(mpi_comm_world,1,ierr)
                end if
            end do
            close(unit=CovIOUnit, iostat=ierr)
            if (ierr/=0) then
                write(*,*) 'Error initializing covariance matrix from file, closing file: ', ierr
                write(*,*) 'I will stop'
                call mpi_abort(mpi_comm_world,1,ierr)
            end if
        end if
        call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
    else
        if (EnsembleRank>NotWorkingMember) then
            BaseIndex=EnsembleRank
        else
            BaseIndex=EnsembleRank+1
        end if
        
        DO jn=1, jptra  ! global loop on tracers to read restart
            write(FileNameBase,'(A31,I3.3,A25)') DirName//'/BASE', BaseIndex, '.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc'
            CALL readnc_slice_double(FileNameBase, 'TRN'//trim(ctrcnm(jn)), BaseMember(:,:,:,jn) )
        end do
        
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
    end if
    
END SUBROUTINE
