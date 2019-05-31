subroutine WriteCov1Seik(FileNameCov)
    use myalloc
    use mpi
    
    implict none
    CHARACTER(LEN=53), intent (in) :: FileNameCov
    integer :: ierr, jn
    
    open(unit=UnitSEIK, file=FileNameCov, form='formatted', iostat=ierr, action='write', access='sequential',status='replace')
    if (ierr/=0) then
        write(*,*) 'Error opening file for writing covariance matrix: ', ierr
        write(*,*) 'I will stop'
        call mpi_abort(mpi_comm_world,1,ierr)
    end if
    do jn=1,SeikDim
        write(UnitSEIK,*,iostat=ierr) CovSeik1(:,jn)
        if (ierr/=0) then
            write(*,*) 'Error writing covariance matrix, while writing line', jn, 'error:', ierr
            write(*,*) 'I will stop'
            call mpi_abort(mpi_comm_world,1,ierr)
        end if
    end do
    close(unit=UnitSEIK, iostat=ierr)
    if (ierr/=0) then
        write(*,*) 'Error closing file of covariance matrix: ', ierr
        write(*,*) 'I will stop'
        call mpi_abort(mpi_comm_world,1,ierr)
    end if 

end subroutine
