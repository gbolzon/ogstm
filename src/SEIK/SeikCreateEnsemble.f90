subroutine SeikCreateEnsemble()
    Use myalloc
    use mpi
    
    implicit none
    integer :: ierr, indexi
    
    if (UseHighOrder) then
        if (MyRank==0) then
            if (EnsembleRank==NotWorkingMember) then
                call dpotrf( 'U', SeikDim, CovSeik1, SeikDim, ierr )
                if (ierr.ne.0) error stop 'Choleky failed!'
                
                ChangeBaseSeik=0.0d0
                do indexi=1,SeikDim 
                    if (indexi>=NotWorkingMember) then
                        ChangeBaseSeik(i,i+1)=1.0d0
                    else
                        ChangeBaseSeik(i,i)=1.0d0
                    end if
                end do
                call dtrtrs( 'U', 'N', 'N', SeikDim, SeikDim+1, CovSeik1, SeikDim, ChangeCoefSeik, SeikDim, ierr)
                if (ierr.ne.0) error stop 'Sampling inversion failed'
            
                call Sampling(CovSeik1, SeikDim, ChangeBaseSeik, ierr)
                call MPI_Scatter(ChangeBaseSeik, SeikDim, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
            else
            
            
                call MPI_Scatter(0, 0, mpi_real8, ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
            end if
        end if
        
        if (.not.(EnsembleRank==NotWorkingMember)) then
            
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
    trnEnsembleWeighted=trnEnsembleWeighted*Weight
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
