subroutine PCAReadMatrices
    use myalloc
    use StringSEIK

    implicit none
    
    integer :: indexi
    logical :: ExistSVD
        
    call LoadMatrix(HistoryForVar,nHistoryForVar,"REDUCED_BASE/PCA/SAVES/HistVar-"//int2str(MyRank,4)//".dat")
    
    INQUIRE(FILE="REDUCED_BASE/PCA/SAVES/HistSVD-"//int2str(MyRank,4)//".dat", EXIST=ExistSVD)
    
    if (ExistSVD) then
        call LoadMatrix(HistoryForSVD,nHistoryForSVD,"REDUCED_BASE/PCA/SAVES/HistSVD-"//int2str(MyRank,4)//".dat")
    else
        do indexi=0, nHistoryForSVD/nHistoryForSVDpart-1
            call LoadMatrix(HistoryForSVD(:,:,:,:,(indexi*nHistoryForSVDpart+1):((indexi+1)*nHistoryForSVDpart)), & 
                nHistoryForSVDpart,"REDUCED_BASE/PCA/SAVES/HistSVDpart"//int2str(indexi,2)//"-"//int2str(MyRank,4)//".dat")
        end do
        
        call SaveMatrix(HistoryForSVD,nHistoryForSVD,"REDUCED_BASE/PCA/SAVES/HistSVD-"//int2str(MyRank,4)//".dat")
    end if
end subroutine    

subroutine LoadMatrix(LoadingMatrix, nColumns, FileName)
    use myalloc
    use mpi

    implicit none

    integer, intent(in) :: nColumns
    double precision, dimension(jpk,jpj,jpi,jptra,nColumns), intent(out) :: LoadingMatrix
    character(len=*), intent(in) :: FileName
    integer :: ierr

    open(unit=UnitSEIK, file=FileName, form='unformatted', iostat=ierr, &
        action='read', access='stream',status='old')
    if (ierr/=0) then
        write(*,*) 'Error opening file while loading matrix: ', ierr
        write(*,*) 'File: ', FileName
        write(*,*) 'I will stop'
        call mpi_abort(mpi_comm_world,1,ierr)
    end if
    
    read(UnitSEIK,iostat=ierr) LoadingMatrix
    if (ierr/=0) then
        write(*,*) 'Error reading file while loading matrix: ', ierr
        write(*,*) 'File: ', FileName
        write(*,*) 'I will stop'
        call mpi_abort(mpi_comm_world,1,ierr)
    end if
    
    close(UnitSEIK,iostat=ierr)
    if (ierr/=0) then
        write(*,*) 'Error closing file while loading matrix: ', ierr
        write(*,*) 'File: ', FileName
        !write(*,*) 'I will stop'
        !call mpi_abort(mpi_comm_world,1,ierr)
    end if
end subroutine
        

 
