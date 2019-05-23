subroutine WriteBaseSeik(DateString)
    use myalloc
    !use mpi
    
    implicit none
    
    character(len=17), intent(in) :: DateString
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
            FileNameCov=DirName//'/COV1.'//DateString//'.csv'
            
            open(unit=UnitSEIK, file=FileNameCov, form='formatted', iostat=ierr, action='write', access='sequential',status='replace')
            if (ierr/=0) then
                write(*,*) 'Error opening file for writing covariance matrix: ', ierr
                !write(*,*) 'I will stop'
                !call mpi_abort(mpi_comm_world,1,ierr)
            end if
            do jn=1,SeikDim
                write(UnitSEIK,*,iostat=ierr) CovSeik1(:,jn)
                if (ierr/=0) then
                    write(*,*) 'Error writing covariance matrix, while writing line', jn, 'error:', ierr
                    !write(*,*) 'I will stop'
                    !call mpi_abort(mpi_comm_world,1,ierr)
                end if
            end do
            close(unit=UnitSEIK, iostat=ierr)
            if (ierr/=0) then
                write(*,*) 'Error closing file of covariance matrix: ', ierr
                !write(*,*) 'I will stop'
                !call mpi_abort(mpi_comm_world,1,ierr)
            end if
        end if
    else
        if (EnsembleRank>NotWorkingMember) then
            BaseIndex=EnsembleRank
        else
            BaseIndex=EnsembleRank+1
        end if
        
        trcwriSeik(DateString, BaseIndex, DirName//'/')
    end if
    
end program
