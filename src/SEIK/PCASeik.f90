subroutine PCASeik
    use mpi
    use myalloc
    use StatSEIK
    
    implicit none
    
    double precision, parameter :: Threshold=1.0d-6
    
    double precision, dimension(nHistoryForSVD,SpaceDim) :: PCAMatrix
    double precision, dimension(nHistoryForSVD) :: tempvec, workvec
    double precision, dimension(nHistoryForSVD,SpaceDim) :: LogMatrix
    double precision , dimension(SpaceDim) :: PCAVar
    double precision , dimension(SpaceDim) :: PCASTD
    double precision, dimension(nHistoryForSVD,nHistoryForSVD) :: LogMatrix2, LogMatrix2part
    integer :: indexi, ierr
    
    double precision, dimension(nHistoryForSVD*nHistoryForSVD) :: work, iwork
    integer dimension(200) :: isuppz
    
    double precision, dimension(nHistoryForSVD) :: eigenvalues
    double precision, dimension(nHistoryForSVD,100) :: eigenvectors
    integer :: neigenvalues
    
    PCAMatrix=transpose(reshape(HistoryForSVD,(/SpaceDim,nHistoryForSVD/)))
    !PCAVar=CalcVar(PCAMatrix, nHistoryForSVD, SpaceDim, CalcMean(PCAMatrix, nHistoryForSVD, SpaceDim))
    
    LogMatrix=0.0d0
    PCASTD=0.0d0
    do indexi=1, SpaceDim
        PCAVar(indexi)=CalcVar(PCAMatrix(:,indexi), nHistoryForSVD, CalcMean(PCAMatrix(:,indexi), nHistoryForSVD), workvec)
        if (PCAVar(indexi)>Threshold) then
            tempvec=log(PCAMatrix(:,indexi))
            tempvec=tempvec-CalcMean(tempvec, nHistoryForSVD)
            PCAVar(indexi)=CalcVar(tempvec, nHistoryForSVD,0.0d0,workvec)
            if (PCAVar(indexi)>1.0d0/ModelErrorDiag1(indexi)) then
                PCASTD=sqrt(PCAVar(indexi)
                LogMatrix(:,indexi)=tempvec/PCASTD)
            end if
        end if
    end do
    
    LogMatrix2part=matmul(LogMatrix,transpose(LogMatrix))
    call mpi_reduce(LogMatrix2part,LogMatrix2, nHistoryForSVD*nHistoryForSVD, mpi_real8, mpi_sum, 0, LocalComm, ierr)
    
    if (MyRank==0) then
        call dsyevr("V", "I", "U", nHistoryForSVD, LogMatrix2, nHistoryForSVD, 0.0d0, 0.0d0, nHistoryForSVD-99, nHistoryForSVD, & 
            1.0d4, neigenvalues, eigenvalues, eigenvectors, nHistoryForSVD, &
            isuppz, work, nHistoryForSVD*nHistoryForSVD, iwork, nHistoryForSVD*nHistoryForSVD, ierr)
        if (ierr/=0) then
            write(*,*) "something wrong with PCA. I will stop"
            call mpi_abort(mpi_comm_world,1,ierr)
        end if
        
        continuare....
        ricorda di rimoltiplicare per PCASTD
    
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

