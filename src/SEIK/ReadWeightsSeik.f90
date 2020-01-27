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
            open(unit=UnitSEIK, file=FileNameCov, form='formatted', iostat=ierr, action='read', access='sequential',status='old')
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

            read(UnitSEIK,*,iostat=ierr) HighOrderMatrix !verifica la costruzione di TTT. davvero si divide per i pesi???
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
            CALL readnc_slice_double(FileNameBase, 'TRN'//trim(ctrcnm(jn)), BaseMember(:,:,:,jn) ) ! This routine should be rewritten...
            BaseMember(:,:,:,jn) = BaseMember(:,:,:,jn) * tmask !*SeikMask
        end do
        
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
    end if
    
END SUBROUTINE
