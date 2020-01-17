module StatSEIK
implicit none

contains
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


subroutine PCASeik !ATTENTION! This routine need a correction: It is necessary to put ghost cells to zero!
    use mpi
    use myalloc
    use StatSEIK
    
    implicit none
    
    double precision, parameter :: Threshold=1.0d-8
    integer, parameter :: maxNeigenvectors=100
    
    !double precision, dimension(nHistoryForSVD,SpaceDim) :: PCAMatrix
    !double precision, dimension(SpaceDim,nHistoryForSVD) :: PCAMatrixT
    double precision, dimension(nHistoryForSVD) :: tempvec, workvec
    !double precision, dimension(nHistoryForSVD,SpaceDim) :: LogMatrix, NormalizedLogMatrix
    double precision , dimension(SpaceDim) :: PCAVar
    !double precision , dimension(SpaceDim) :: PCASTD
    double precision, dimension(nHistoryForSVD,nHistoryForSVD) :: NormalizedLogMatrix2, NormalizedLogMatrix2part
    integer :: indexi, indexj, indexk, indexc, indexn, ierr
    
    double precision, dimension(nHistoryForSVD*nHistoryForSVD) :: work, iwork
    integer, dimension(2*maxNeigenvectors) :: isuppz
    double precision dlamch
    
    double precision, dimension(nHistoryForSVD) :: eigenvalues
    double precision, dimension(nHistoryForSVD,maxNeigenvectors) :: eigenvectors
    integer :: neigenvalues
    
    double precision, dimension(:,:), allocatable :: Transformation
    double precision :: temptime

    if(myrank==0) write(*,*) "Starting PCA"
    temptime=mpi_wtime()
    !call mpi_barrier(EnsembleComm, ierr) qui in teoria abbiamo un solo member, che senso ha chiamare un barrier tra i members?
    if(myrank==0) write(*,*) "building matrix"

if (.false.) then !true to calcolate mean, false for svd
    !where (HistoryForSVD>1.0d-8)
    !    HistoryForSVD=log(HistoryForSVD)
    !elsewhere
    !    HistoryForSVD=log(1.0d-8)
    !end where
    do indexn=1, nHistoryForSVD
        do indexc=1, jptra
            do indexi=1, jpi
                do indexj=1, jpj
                    do indexk=1, jpk
                        if (HistoryForSVD(indexk,indexj,indexi,indexc,indexn)>1.0d-11) then
                            HistoryForSVD(indexk,indexj,indexi,indexc,indexn)=log(HistoryForSVD(indexk,indexj,indexi,indexc,indexn))
                        else
                            HistoryForSVD(indexk,indexj,indexi,indexc,indexn)=log(1.0d-11)
                        end if
                    end do
                end do
            end do
        end do
    end do
    call mpi_barrier(EnsembleComm, ierr) !Errore, questi sono tutti LocalComm???
    if (myrank==0) write(*,*) "log"
    trn=sum(HistoryForSVD,5)
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "sum"
    trn=trn/nHistoryForSVD
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "mean"
    BaseMember=exp(trn)    
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "exp"
    !call trcwri("20000101-00:00:00")
    call trcwriSeik("20000101-00:00:00",-1, 'REDUCED_BASE/PCA/', BaseMember)

else    
    
    write(*,*) "r=", myrank, "sum SeikMask=", sum(SeikMask)
    do indexi=1, nHistoryForSVD
        do indexj=1, jptra
            where (SeikMask==0) HistoryForSVD(:,:,:,indexj, indexi)=0.0d0 !instead of bfmmask
            !HistoryForSVD(:,:,:,indexj, indexi)=HistoryForSVD(:,:,:,indexj, indexi)*bfmmask
        end do
        where (SeikTrcMask==0) HistoryForSVD(:,:,:,:, indexi)=0.0d0
        call CutCellsTracer(HistoryForSVD(:,:,:,:, indexi))
    end do
    PCAVar=0.0d0
    
    PCAMatrixT=reshape(HistoryForSVD,(/SpaceDim,nHistoryForSVD/))
    PCAMatrix=transpose(PCAMatrixT)
    !PCAVar=CalcVar(PCAMatrix, nHistoryForSVD, SpaceDim, CalcMean(PCAMatrix, nHistoryForSVD, SpaceDim))
    
    !LogMatrix=0.0d0
    !PCASTD=0.0d0
    !where (PCAMatrix<1.0d-6) PCAMatrix=1.0d-6

    if (.true.) then ! checking
        do indexi=1, SpaceDim
            PCAVar(indexi)=CalcVar(PCAMatrix(:,indexi), nHistoryForSVD, CalcMean(PCAMatrix(:,indexi), nHistoryForSVD), workvec)
        end do
        BaseMember=reshape(PCAVar,(/jpk,jpj,jpi,jptra/))
        !where (BaseMember.le.Threshold) BaseMember=1.0d20
        call trcwriSeik('PCAVar78901234567', -1, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)
        
        BaseMember=sqrt(BaseMember)
        !where (BaseMember.le.Threshold) BaseMember=1.0d20
        call trcwriSeik('PCAStd78901234567', -1, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)

        !do indexi=1, jptra
        !    BaseMember(:,:,:,indexi)=bfmmask
	!end do
        !call trcwriSeik('bfmmask8901234567', 101, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)

        !BaseMember=reshape(HistoryForSVD(:,:,:,:, 10),(/jpk,jpj,jpi,jptra/))
        !call trcwriSeik('Hist5678901234567', 101, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)

        call mpi_barrier(mpi_comm_World,ierr)
        if (myrank==0) write(*,*) "checks written"
        !call mpi_abort(mpi_comm_world,1,ierr)
    end if
    
    where (PCAMatrix<CutOffValue) PCAMatrix=CutOffValue
    
    do indexi=1, SpaceDim
        !PCAVar(indexi)=CalcVar(PCAMatrix(:,indexi), nHistoryForSVD, CalcMean(PCAMatrix(:,indexi), nHistoryForSVD), workvec)
        !if (PCAVar(indexi)>Threshold) then
            tempvec=log(PCAMatrix(:,indexi))
            tempvec=tempvec-CalcMean(tempvec, nHistoryForSVD)
            PCAVar(indexi)=CalcVar(tempvec, nHistoryForSVD,0.0d0,workvec)
            !if (PCAVar(indexi)>1.0d0/ModelErrorDiag1(indexi)) then
            if (PCAVar(indexi)>1.0d-8) then
                !PCASTD(indexi)=sqrt(PCAVar(indexi))
                PCAVar(indexi)=sqrt(PCAVar(indexi))
                !LogMatrix(:,indexi)=tempvec
                !NormalizedLogMatrix(:,indexi)=tempvec/PCASTD(indexi)
                !NormalizedLogMatrix(:,indexi)=tempvec/sqrt(PCAVar(indexi))
                PCAMatrix(:,indexi)=tempvec/PCAVar(indexi)
            else
                PCAMatrix(:,indexi)=0.0d0
                PCAVar(indexi)=0.0d0
            end if
        !else
            !PCAMatrix(:,indexi)=0.0d0
            !PCAVar(indexi)=0.0d0
        !end if
    end do
    
    if (.true.) then
        BaseMember=reshape(PCAVar,(/jpk,jpj,jpi,jptra/))
        call trcwriSeik('PCALogStd01234567', -1, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)
        call mpi_barrier(mpi_comm_World,ierr)
    end if

    if(myrank==0) write(*,*) "reducing matrix"    

    !NormalizedLogMatrix2part=matmul(NormalizedLogMatrix,transpose(NormalizedLogMatrix))
    PCAMatrixT=transpose(PCAMatrix)
    NormalizedLogMatrix2part=matmul(PCAMatrix,PCAMatrixT)
    call mpi_reduce(NormalizedLogMatrix2part,NormalizedLogMatrix2, nHistoryForSVD*nHistoryForSVD, mpi_real8, mpi_sum, 0, LocalComm, ierr)
    
    if(myrank==0) write(*,*) "svd"

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
    
    if(myrank==0) write(*,*) "preparing and writing base"

    do indexi=1, SpaceDim
        if (UseMaxVarSEIK) then
            if (PCAVar(indexi)>sqrt(MaxVarSEIK)) PCAVar(indexi)=sqrt(MaxVarSEIK)
        end if
        PCAMatrix(:, indexi)=PCAMatrix(:, indexi)*PCAVar(indexi)
    end do

    do indexi=1, neigenvalues
    
        !reusing pcavar instead of defining a new temp array
        PCAVar=matmul(Transformation(:,indexi),PCAMatrix)
        BaseMember=reshape(PCAVar,(/jpk,jpj,jpi,jptra/))
        BaseMember=BaseMember/sqrt(dble(nHistoryForSVD))
        call trcwriSeik('20121231-00:00:00', indexi, 'REDUCED_BASE/PCA/', BaseMember)
        
    end do
    
    deallocate(Transformation)

end if
    
    if(myrank==0) write(*,*) "End of PCA in ", mpi_wtime()-temptime, " seconds"

end subroutine

subroutine PCASeikOld !ATTENTION! This routine need a correction: It is necessary to put ghost cells to zero!
    use mpi
    use myalloc
    use StatSEIK
    
    implicit none
    
    double precision, parameter :: Threshold=1.0d-8
    integer, parameter :: maxNeigenvectors=100
    
    !double precision, dimension(nHistoryForSVD,SpaceDim) :: PCAMatrix
    !double precision, dimension(SpaceDim,nHistoryForSVD) :: PCAMatrixT
    double precision, dimension(nHistoryForSVD) :: tempvec, workvec
    !double precision, dimension(nHistoryForSVD,SpaceDim) :: LogMatrix, NormalizedLogMatrix
    double precision , dimension(SpaceDim) :: PCAVar
    !double precision , dimension(SpaceDim) :: PCASTD
    double precision, dimension(nHistoryForSVD,nHistoryForSVD) :: NormalizedLogMatrix2, NormalizedLogMatrix2part
    integer :: indexi, indexj, ierr
    
    double precision, dimension(nHistoryForSVD*nHistoryForSVD) :: work, iwork
    integer, dimension(2*maxNeigenvectors) :: isuppz
    double precision dlamch
    
    double precision, dimension(nHistoryForSVD) :: eigenvalues
    double precision, dimension(nHistoryForSVD,maxNeigenvectors) :: eigenvectors
    integer :: neigenvalues
    
    double precision, dimension(:,:), allocatable :: Transformation
    double precision :: temptime

    if(myrank==0) write(*,*) "Starting PCA"
    temptime=mpi_wtime()
    if(myrank==0) write(*,*) "building matrix"

if (.false.) then !true to calcolate mean, false for svd
    where (HistoryForSVD>1.0d-8)
        HistoryForSVD=log(HistoryForSVD)
    elsewhere
        HistoryForSVD=log(1.0d-8)
    end where
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "log"
    trn=sum(HistoryForSVD,4)
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "sum"
    trn=trn/nHistoryForSVD
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "mean"
    BaseMember=exp(trn)    
    call mpi_barrier(EnsembleComm, ierr)
    if (myrank==0) write(*,*) "exp"
    !call trcwri("20000101-00:00:00")
    call trcwriSeik("20000101-00:00:00",-1, 'REDUCED_BASE/PCA/', BaseMember)

else    

    do indexi=1, nHistoryForSVD
        do indexj=1, jptra
            where (SeikMask==0) HistoryForSVD(:,:,:,indexj, indexi)=0.0d0 !instead of bfmmask
            !HistoryForSVD(:,:,:,indexj, indexi)=HistoryForSVD(:,:,:,indexj, indexi)*bfmmask
        end do
        call CutCellsTracer(HistoryForSVD(:,:,:,:, indexi))
    end do
    PCAVar=0.0d0
    
    PCAMatrixT=reshape(HistoryForSVD,(/SpaceDim,nHistoryForSVD/))
    PCAMatrix=transpose(PCAMatrixT)
    !PCAVar=CalcVar(PCAMatrix, nHistoryForSVD, SpaceDim, CalcMean(PCAMatrix, nHistoryForSVD, SpaceDim))
    
    !LogMatrix=0.0d0
    !PCASTD=0.0d0
    where (PCAMatrix<1.0d-6) PCAMatrix=1.0d-6

    if (.true.) then ! checking
        do indexi=1, SpaceDim
            PCAVar(indexi)=CalcVar(PCAMatrix(:,indexi), nHistoryForSVD, CalcMean(PCAMatrix(:,indexi), nHistoryForSVD), workvec)
        end do
        BaseMember=reshape(PCAVar,(/jpk,jpj,jpi,jptra/))
        where (BaseMember.le.Threshold) BaseMember=1.0d20
        call trcwriSeik('PCVar678901234567', 101, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)

        !do indexi=1, jptra
        !    BaseMember(:,:,:,indexi)=bfmmask
	!end do
        !call trcwriSeik('bfmmask8901234567', 101, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)

        !BaseMember=reshape(HistoryForSVD(:,:,:,:, 10),(/jpk,jpj,jpi,jptra/))
        !call trcwriSeik('Hist5678901234567', 101, 'REDUCED_BASE/PCA/EXTRA/',BaseMember)

        call mpi_barrier(mpi_comm_World,ierr)
        !call mpi_abort(mpi_comm_world,1,ierr)
    end if
    
    do indexi=1, SpaceDim
        PCAVar(indexi)=CalcVar(PCAMatrix(:,indexi), nHistoryForSVD, CalcMean(PCAMatrix(:,indexi), nHistoryForSVD), workvec)
        if (PCAVar(indexi)>Threshold) then
            tempvec=log(PCAMatrix(:,indexi))
            tempvec=tempvec-CalcMean(tempvec, nHistoryForSVD)
            PCAVar(indexi)=CalcVar(tempvec, nHistoryForSVD,0.0d0,workvec)
            if (PCAVar(indexi)>1.0d0/ModelErrorDiag1(indexi)) then
                !PCASTD(indexi)=sqrt(PCAVar(indexi))
                PCAVar(indexi)=sqrt(PCAVar(indexi))
                !LogMatrix(:,indexi)=tempvec
                !NormalizedLogMatrix(:,indexi)=tempvec/PCASTD(indexi)
                !NormalizedLogMatrix(:,indexi)=tempvec/sqrt(PCAVar(indexi))
                PCAMatrix(:,indexi)=tempvec/PCAVar(indexi)
            else
                PCAMatrix(:,indexi)=0.0d0
                PCAVar(indexi)=0.0d0
            end if
        else
            PCAMatrix(:,indexi)=0.0d0
            PCAVar(indexi)=0.0d0
        end if
    end do

    if(myrank==0) write(*,*) "reducing matrix"    

    !NormalizedLogMatrix2part=matmul(NormalizedLogMatrix,transpose(NormalizedLogMatrix))
    PCAMatrixT=transpose(PCAMatrix)
    NormalizedLogMatrix2part=matmul(PCAMatrix,PCAMatrixT)
    call mpi_reduce(NormalizedLogMatrix2part,NormalizedLogMatrix2, nHistoryForSVD*nHistoryForSVD, mpi_real8, mpi_sum, 0, LocalComm, ierr)
    
    if(myrank==0) write(*,*) "svd"

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
    
    if(myrank==0) write(*,*) "preparing and writing base"

    do indexi=1, SpaceDim
        PCAMatrix(:, indexi)=PCAMatrix(:, indexi)*PCAVar(indexi)
    end do

    do indexi=1, neigenvalues
    
        !reusing pcavar instead of defining a new temp array
        PCAVar=matmul(Transformation(:,indexi),PCAMatrix)
        BaseMember=reshape(PCAVar,(/jpk,jpj,jpi,jptra/))
        BaseMember=BaseMember/sqrt(dble(nHistoryForSVD))
        call trcwriSeik('19990101-00:00:00', indexi, 'REDUCED_BASE/PCA/', BaseMember)
        
    end do
    
    deallocate(Transformation)

end if
    
    if(myrank==0) write(*,*) "End of PCA in ", mpi_wtime()-temptime, " seconds"

end subroutine

