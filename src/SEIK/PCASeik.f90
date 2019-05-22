subroutine PCASeik
    use mpi
    use myalloc
    use StatSEIK
    
    implicit none
    
    double precision, parameter :: Threshold=1.0d-8
    integer, parameter :: maxNeigenvectors
    
    double precision, dimension(nHistoryForSVD,SpaceDim) :: PCAMatrix
    double precision, dimension(nHistoryForSVD) :: tempvec, workvec
    double precision, dimension(nHistoryForSVD,SpaceDim) :: LogMatrix, NormalizedLogMatrix
    double precision , dimension(SpaceDim) :: PCAVar
    !double precision , dimension(SpaceDim) :: PCASTD
    double precision, dimension(nHistoryForSVD,nHistoryForSVD) :: NormalizedLogMatrix2, NormalizedLogMatrix2part
    integer :: indexi, ierr
    
    double precision, dimension(nHistoryForSVD*nHistoryForSVD) :: work, iwork
    integer dimension(200) :: isuppz
    double precision dlamch
    
    double precision, dimension(nHistoryForSVD) :: eigenvalues
    double precision, dimension(nHistoryForSVD,maxNeigenvectors) :: eigenvectors
    integer :: neigenvalues
    
    double precision, dimension(:,:), allocatable :: Transformation
    
    PCAMatrix=transpose(reshape(HistoryForSVD,(/SpaceDim,nHistoryForSVD/)))
    !PCAVar=CalcVar(PCAMatrix, nHistoryForSVD, SpaceDim, CalcMean(PCAMatrix, nHistoryForSVD, SpaceDim))
    
    LogMatrix=0.0d0
    !PCASTD=0.0d0
    do indexi=1, SpaceDim
        PCAVar(indexi)=CalcVar(PCAMatrix(:,indexi), nHistoryForSVD, CalcMean(PCAMatrix(:,indexi), nHistoryForSVD), workvec)
        if (PCAVar(indexi)>Threshold) then
            tempvec=log(PCAMatrix(:,indexi))
            tempvec=tempvec-CalcMean(tempvec, nHistoryForSVD)
            PCAVar(indexi)=CalcVar(tempvec, nHistoryForSVD,0.0d0,workvec)
            if (PCAVar(indexi)>1.0d0/ModelErrorDiag1(indexi)) then
                !PCASTD(indexi)=sqrt(PCAVar(indexi))
                LogMatrix(:,indexi)=tempvec
                !NormalizedLogMatrix(:,indexi)=tempvec/PCASTD(indexi)
                NormalizedLogMatrix(:,indexi)=tempvec/sqrt(PCAVar(indexi))
            end if
        end if
    end do
    
    NormalizedLogMatrix2part=matmul(NormalizedLogMatrix,transpose(NormalizedLogMatrix))
    call mpi_reduce(NormalizedLogMatrix2part,NormalizedLogMatrix2, nHistoryForSVD*nHistoryForSVD, mpi_real8, mpi_sum, 0, LocalComm, ierr)
    
    if (MyRank==0) then
        call dsyevr("V", "I", "U", nHistoryForSVD, NormalizedLogMatrix2, nHistoryForSVD, 0.0d0, 0.0d0, nHistoryForSVD-maxNeigenvectors+1, nHistoryForSVD, & 
            dlamch('S'), neigenvalues, eigenvalues, eigenvectors, nHistoryForSVD, &
            isuppz, work, nHistoryForSVD*nHistoryForSVD, iwork, nHistoryForSVD*nHistoryForSVD, ierr)
        if (ierr/=0) then
            write(*,*) "something wrong with PCA. I will stop"
            call mpi_abort(mpi_comm_world,1,ierr)
        end if
        
        if (maxNeigenvectors/=neigenvalues) then
            write(*,*) "something strange in the number of eigenvalues!"
            write(*,*) "maxNeigenvectors=", maxNeigenvectors, " neigenvalues=", neigenvalues
        end if
        
        open(unit=UnitSEIK, file='REDUCED_BASE/PCA/eigenvalues.csv', form='formatted', iostat=ierr, action='write', access='sequential',status='replace')
        if (ierr/=0) then
            write(*,*) 'Error opening file for eigenevalues: ', ierr
            write(*,*) 'I will stop'
            call mpi_abort(mpi_comm_world,1,ierr)
        end if
        write(UnitSEIK,*,iostat=ierr) eigenvalues
        if (ierr/=0) then
            write(*,*) 'Error writing eigenvalues:', ierr
            write(*,*) 'I will stop'
            call mpi_abort(mpi_comm_world,1,ierr)
        end if
        close(unit=UnitSEIK, iostat=ierr)
        if (ierr/=0) then
            write(*,*) 'Error closing eigenvalues file: ', ierr
            write(*,*) 'I will stop'
            call mpi_abort(mpi_comm_world,1,ierr)
        end if
        
        do indexi=1, maxNeigenvectors
            if (eigenvalues(indexi)>1.0d-12) then
                exit
            else
                neigenvalues=neigenvalues-1
            end if
        end do
        
        call mpi_bcast(neigenvalues,1, MPI_int, 0, LocalComm, ierr)
        
        allocate(Transformation(nHistoryForSVD,neigenvalues))
        Transformation(:,:)=eigenvectors(:,maxNeigenvectors:maxNeigenvectors-neigenvalues+1:-1)
        
        call mpi_bcast(Transformation,nHistoryForSVD*neigenvalues, mpi_real8, 0, LocalComm, ierr)
    
    else
        
        call mpi_bcast(neigenvalues,1, MPI_int, 0, LocalComm, ierr)
        
        allocate(Transformation(nHistoryForSVD,neigenvalues))
        
        call mpi_bcast(Transformation,nHistoryForSVD*neigenvalues, mpi_real8, 0, LocalComm, ierr)
    
    end if
    
    do indexi=1, neigenvalues
    
        !reusing pcavar instead of defining a new temp array
        PCAVar=matmul(Transformation(:,indexi),LogMatrix)
        BaseMember=reshape(PCAVar,(/jpk,jpj,jpi,jptra/))
        
        call trcwriSeik('19990101-00:00:00', indexi, 'REDUCED_BASE/PCA/')
        
    end do

end subroutine


module StatSEIK
implicit none

conteins
function CalcVar(array, nside, meanval, tempvec)
    implicit none
    
    integer, intent(in) :: nside
    double precision, dimension(nside), intent(in) :: array
    double precision, intent(in) :: meanval
    double precision :: CalcVar
    double precision, dimension(nside), intent(inout) :: tempvec
    
    tempvec=array-meanval
    tempvec=tempvec*tempvec
    CalcVar=sum(tempvec)/nside

end function

function CalcMean(array, nside)
    implicit none
    
    integer, intent(in) :: nside
    double precision, dimension(nside), intent(in) :: array
    double precision:: CalcMean
    
    CalcMean=sum(array)/nside
end function

end module

