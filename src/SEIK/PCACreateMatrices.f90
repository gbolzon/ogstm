subroutine PCACreateMatrices(DateString)
    use myalloc
    use StringSEIK

    implicit none
    
    character(len=17), intent(in) :: DateString
    integer :: year, monthday
        
    if ((PCAFullYear).and.(DATEstring(1:4).eq."2000").and.(DATEstring(10:17).eq."00:00:00")) then
        CounterForVar=CounterForVar+1
        HistoryForVar(:,:,:,:,CounterForVar)=trn
        if (CounterForVar==nHistoryForVar) then
            call SaveMatrix(HistoryForVar,nHistoryForVar,"REDUCED_BASE/PCA/SAVES/HistVar-"//int2str(MyRank,4)//".dat")
        end if
    end if
    
    if (DATEstring(10:17).eq."00:00:00") then
        year=str2int(DATEstring(1:4))
        monthday=str2int(DATEstring(5:8))
        if (((year.ge.2014).and.(year.le.2009).and.(monthday.le.116)) .or. &
            ((year.ge.2013).and.(year.le.2008).and.(monthday.ge.1217))) then
            
            CounterForSVDpart=CounterForSVDpart+1
            HistoryForSVDpart(:,:,:,:,CounterForSVDpart)=trn
            if (CounterForSVDpart==nHistoryForSVDpart) then
                CounterForSVDpart=0
                call SaveMatrix(HistoryForSVDpart,nHistoryForSVDpart,"REDUCED_BASE/PCA/SAVES/HistSVDpart"//int2str(SVDpartID,2)//"-"//int2str(MyRank,4)//".dat")
                SVDpartID=SVDpartID+1
            end if
            
            CounterForSVD=CounterForSVD+1
            HistoryForSVD(:,:,:,:,CounterForSVD)=trn
            if (CounterForSVD==nHistoryForSVD) then
                call SaveMatrix(HistoryForSVD,nHistoryForSVD,"REDUCED_BASE/PCA/SAVES/HistSVD-"//int2str(MyRank,4)//".dat")
            end if
        end if
    end if
end subroutine    

subroutine SaveMatrix(SavingMatrix, nColumns, FileName)
    use myalloc
    !use mpi

    implicit none

    integer, intent(in) :: nColumns
    double precision, dimension(jpk,jpj,jpi,jptra,nColumns), intent(in) :: SavingMatrix
    character(len=*), intent(in) :: FileName
    integer :: ierr

    open(unit=UnitSEIK, file=FileName, form='unformatted', iostat=ierr, &
        action='write', access='stream',status='new')
    if (ierr/=0) then
        write(*,*) 'Error opening file while saving matrix: ', ierr
        write(*,*) 'File: ', FileName
        !write(*,*) 'I will stop'
        !call mpi_abort(mpi_comm_world,1,ierr)
    end if
    
    write(UnitSEIK,iostat=ierr) SavingMatrix
    if (ierr/=0) then
        write(*,*) 'Error writing file while saving matrix: ', ierr
        write(*,*) 'File: ', FileName
        !write(*,*) 'I will stop'
        !call mpi_abort(mpi_comm_world,1,ierr)
    end if
    
    close(UnitSEIK,iostat=ierr)
    if (ierr/=0) then
        write(*,*) 'Error closing file while saving matrix: ', ierr
        write(*,*) 'File: ', FileName
        !write(*,*) 'I will stop'
        !call mpi_abort(mpi_comm_world,1,ierr)
    end if
end subroutine
        

