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
