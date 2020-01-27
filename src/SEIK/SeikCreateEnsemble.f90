subroutine SeikCreateEnsemble()
    Use myalloc
    use mpi
    
    implicit none
    integer :: ierr, indexi, neigenvalues
    double precision dlamch
    
    if (UseHighOrder) then

if (.false.) then
if (MyRank==0) then
    if (EnsembleRank==NotWorkingMember) then
        call dpotrf( 'U', SeikDim, CovSeik1, SeikDim, ierr )
        if (ierr.ne.0) error stop 'Cholesky failed!'
        
        ChangeBaseSeik=0.0d0
        do indexi=1,SeikDim 
            if (indexi>=NotWorkingMember) then
                ChangeBaseSeik(indexi,indexi+1)=1.0d0
            else
                ChangeBaseSeik(indexi,indexi)=1.0d0
            end if
        end do
        call dtrtrs( 'U', 'N', 'N', SeikDim, SeikDim+1, CovSeik1, SeikDim, ChangeBaseSeik, SeikDim, ierr)
        if (ierr.ne.0) error stop 'Sampling inversion failed'
    
        call MPI_Scatter(ChangeBaseSeik, SeikDim, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
    else
        call MPI_Scatter(0, 0, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        
        call MPI_Bcast(ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
        TempVecSeik=matmul(Lseik,ChangeCoefSeik)
        BaseMember=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))
        
    end if
end if

if (.not.(EnsembleRank==NotWorkingMember)) then
    
end if
end if
        
        !da parallelizzare bene
        if (EnsembleRank==NotWorkingMember) then
            do indexi=1,SeikDim 
                BaseMember=reshape(LSeik(:,indexi), (/jpk,jpj,jpi,jptra/)
                where (SeikTrcMask==0) BaseMember=0.0d0
                call CutCellsTracer(BaseMember)
                LSeik(:,indexi)=reshape(BaseMember, (/SpaceDim/)
            end do
            
            LSeikT=transpose(LSeik)
            
            if (MyRank==0) then
                if (UseCholesky) then
                    call dpotrf( 'L', SeikDim, CovSeik1, SeikDim, ierr )
                    if (ierr.ne.0) error stop 'Cholesky failed!'
                else
                    call dsyevr("V", "A", "U", SeikDim, CovSeik1, SeikDim, 0.0d0, 0.0d0,0.0d0, 0.0d0, &
                        dlamch('S'), neigenvalues, eigenvalues, eigenvectors, SeikDim, &
                        isuppz, work, SeikDim*SeikDim, iwork, SeikDim*SeikDim, ierr)

                    if (ierr/=0) then
                        write(*,*) "something wrong with svd. I will stop"
                        call mpi_abort(mpi_comm_world,1,ierr)
                    end if

                    if (SeikDim/=neigenvalues) then
                        write(*,*) "something strange in the number of eigenvalues!"
                        write(*,*) "maxNeigenvectors=", maxNeigenvectors, " neigenvalues=", neigenvalues
                    end if
                    
                    do indexi=1, SeikDim
                        CovSeik1(:,indexi)=eigenvectors(:,indexi)/sqrt(eigenvalues(indexi))
                    end do
                    CovSeik1=Inverse(CovSeik1)
                    CovSeik1=MatMul(eigenvectos,CovSeik1)
               end if
            end if
            
            call MPI_Bcast(CovSeik1, SeikDim*SeikDim, mpi_real8, 0, LocalComm, ierr)
            
            if (UseCholesky) then
                call dtrtrs( 'L', 'N', 'N', SeikDim, SpaceDim, CovSeik1, SeikDim, LSeikT, SeikDim, ierr)
                if (ierr.ne.0) error stop 'Sampling inversion failed'
            else
                LSeikT=MatMul(CovSeik1,LSeikT)
            end if
    
            do indexi=1,SpaceDim
                TempSliceSeik=LSeikT(:,indexi)
                TempVecSeik(indexi)=norm2(TempSliceSeik)
                if (TempVecSeik(indexi)>1.0d-4) then
                    LSeikT(:,indexi)=TempSliceSeik/TempVecSeik(indexi)
                else
                    LSeikT(:,indexi)=0.0d0
                    TempVecSeik(indexi)=0.0d0
                end if
            end do
            
            LSeik=transpose(LSeikT)
            CovSeik1=matmul(LSeikT,LSeik)
            call mpi_reduce(CovSeik1,SvdMatrix, SeikDim*SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
        
            if (MyRank==0) then
                call dsyevr("V", "A", "U", SeikDim, SvdMatrix, SeikDim, 0.0d0, 0.0d0,0.0d0, 0.0d0, & 
                    dlamch('S'), neigenvalues, eigenvalues, eigenvectors, SeikDim, &
                    isuppz, work, SeikDim*SeikDim, iwork, SeikDim*SeikDim, ierr)
                    
                if (ierr/=0) then
                    write(*,*) "something wrong with svd. I will stop"
                    call mpi_abort(mpi_comm_world,1,ierr)
                end if
                
                if (SeikDim/=neigenvalues) then
                    write(*,*) "something strange in the number of eigenvalues!"
                    write(*,*) "maxNeigenvectors=", maxNeigenvectors, " neigenvalues=", neigenvalues
                end if
                
                CovSeik1=eigenvectors(:,SeikDim:+1:-1)
                
                call mpi_bcast(CovSeik1,SeikDim*SeikDim, mpi_real8, 0, LocalComm, ierr)
            
            else
            
                call mpi_bcast(CovSeik1,SeikDim*SeikDim, mpi_real8, 0, LocalComm, ierr)
            
            end if
        
            do indexi=1, SpaceDim
                if (UseMaxVarSEIK) then
                    if (TempVecSeik(indexi)>sqrt(MaxVarSEIK)) TempVecSeik(indexi)=sqrt(MaxVarSEIK)
                end if
                LSeikT(:, indexi)=LSeikT(:, indexi)*TempVecSeik(indexi)
            end do

            do indexi=1, SeikDim
                LSeik(:,indexi)=matmul(CovSeik1(:,indexi),LSeikT)              
            end do
        
            call SamplingHighOrder(CovSeik1, SeikDim, ChangeBaseSeik, ierr)
        
        end if
        
        
        
    end if
    
    if (MyRank==0) then
        if (EnsembleRank==NotWorkingMember) then
            call Sampling(CovSeik1, SeikDim, ChangeBaseSeik, ierr)
            call MPI_Scatter(ChangeBaseSeik, SeikDim, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        else
            call MPI_Scatter(0, 0, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
        end if
    end if
    call MPI_Bcast(ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    TempVecSeik=matmul(Lseik,ChangeCoefSeik)
    BaseMember=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))

if (.false.) then
! questa parte va rimossa, perche' inutile durante il main loop, ma non ho tempo ora.
! L'unico momento in cui serve e' alla creazione del primo ensamble (se la PCA non è gia'
! stata fatta con UseMaxVarSEIK). Bisognerebbe metterla su ReadBaseSeik.
! Per il momento la tolgo, tanto al massimo avrò un po' di instabilita' al primo timestep.
if (UseMaxVarSEIK) then
    trnEnsembleWeighted=BaseMember**2
    trnEnsembleWeighted=trnEnsembleWeighted*SeikWeight
    call MPI_AllReduce(trnEnsembleWeighted, trnVariance, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    where (trnVariance>MaxVarSEIK)
        trnEnsemble=sqrt(MaxVarSEIK/trnVariance)
    elsewhere
        trnEnsemble=1.0d0
    end where
    BaseMember=BaseMember*trnEnsemble
end if
end if

    BaseMember=exp(BaseMember)
    trn=trn*BaseMember
    trb=trn
end subroutine
