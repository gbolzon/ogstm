SUBROUTINE ReadWeightsSeik

    USE myalloc
    use mpi

    IMPLICIT NONE

    INTEGER :: jn, BaseIndex, ierr
    CHARACTER(LEN=20) :: DirName
    CHARACTER(LEN=32) :: FileName

    DirName='REDUCED_BASE/WEIGHTS'
    FileName='REDUCED_BASE/WEIGHTS/WGT_006.csv'
    
    write(FileName,'(A25,I3.3,A4)') DirName//'/WGT_', SeikDim,'.csv'
    if (EnsembleRank==NotWorkingMember) then
        if (myrank==0) then
            open(unit=UnitSEIK, file=FileName, form='formatted', iostat=ierr, action='read', access='sequential',status='old')
            if (ierr/=0) then
                write(*,*) 'Error initializing weights matrix from file: ', ierr
                write(*,*) 'I will stop'
                call mpi_abort(mpi_comm_world,1,ierr)
            end if
            
            read(UnitSEIK,*,iostat=ierr) HighOrderDim
            if (ierr/=0) then
                write(*,*) 'Error initializing weights matrix from file, while reading dimension. Error:', ierr
                write(*,*) 'I will stop'
                call mpi_abort(mpi_comm_world,1,ierr)
            end if

            allocate(HighOrderMatrix(SeikDim+1,HighOrderDim+1))

            read(UnitSEIK,*,iostat=ierr) HighOrderMatrix !verifica la costruzione di TTT. davvero si divide per i pesi??? Si', X T (T^T W^-1 T)^-1 T^T W^-1 = X se X ha media nulla.
            if (ierr/=0) then
                write(*,*) 'Error initializing weights matrix from file, while reading lines. Error:', ierr
                write(*,*) 'I will stop'
                call mpi_abort(mpi_comm_world,1,ierr)
            end if

            close(unit=UnitSEIK, iostat=ierr)
            if (ierr/=0) then
                write(*,*) 'Error initializing weights matrix from file, closing file: ', ierr
                write(*,*) 'I will stop'
                call mpi_abort(mpi_comm_world,1,ierr)
            end if

            AllWeightsSqrt=HighOrderMatrix(:,1)
            AllWeights=AllWeightsSqrt**2
            AllWeightsSqrt1=1.0d0/AllWeightsSqrt
            call TTTSeik_builder()
            call MPI_Scatter(AllWeights,SeikDim+1,MPI_Real8,SeikWeight,1,MPI_Real8,NotWorkingMember,EnsembleComm,ierr)
        end if
        
    else
        if (myrank==0) call MPI_Scatter(0,0,MPI_Real8,SeikWeight,1,MPI_Real8,NotWorkingMember,EnsembleComm,ierr)
    end if
    call MPI_Bcast( SeikWeight, 1, MPI_Real8, 0,LocalComm, ierr)

END SUBROUTINE
