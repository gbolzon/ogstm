subroutine WriteBaseSeik(DateString)
    use myalloc
    !use mpi
    
    implicit none
    
    character(len=17), intent(in) :: DateString
    INTEGER :: jn, BaseIndex, ierr
    CHARACTER(LEN=17) :: DirName
    !CHARACTER(LEN=50) :: FileNameBase
    CHARACTER(LEN=44) :: FileNameCov
    
    DirName='REDUCED_BASE/BASE'
    !FileNameBase='REDUCED_BASE/BASE/BASE001.20111231-15:30:00.N1p.nc'
    !FileNameCov='REDUCED_BASE/BASE/COV1.20111231-15:30:00.csv'
    
    !write(DirName,'(A23,I3.3)') 'REDUCED_BASE/DIMENSION_', SeikDim
    if (EnsembleRank==NotWorkingMember) then
        if (myrank==0) then
            FileNameCov=DirName//'/COV1.'//DateString//'.csv'
            call WriteCov1Seik(FileNameCov)
        end if
    else
        if (EnsembleRank>NotWorkingMember) then
            BaseIndex=EnsembleRank
        else
            BaseIndex=EnsembleRank+1
        end if
        
        call trcwriSeik(DateString, BaseIndex, DirName//'/', BaseMember)
    end if
    
end subroutine
