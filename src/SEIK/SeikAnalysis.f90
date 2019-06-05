subroutine SeikAnalysis
    use mpi
    use myalloc
use StringSEIK    

    implicit none
    
    integer :: ierr
integer jn
    
    ObsBaseMember=ComputedObsSeik*Weight
    call MPI_AllReduce(ObsBaseMember, ComputedObsMean, ObsSpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    ObsBaseMember=ComputedObsSeik-ComputedObsMean
    
    if (EnsembleRank==NotWorkingMember) then
        call mpi_allgatherv(0,0,mpi_real8,HLSeik,MpiCountObs,MpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
        TempObsSeik=reshape(MisfitSeik,(/ObsSpaceDim/))
        TempObsSeik=TempObsSeik*ObsErrorDiag1
        TempSliceSeik=matmul(TempObsSeik,HLSeik)
        if (MyRank==0) then
            call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            call mpi_gatherv(0,0,mpi_real8,HLTR1HL,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)

write(*,*) "controllo finale", EnsembleRank, MyRank, CovSeik1
write(*,*) "controllo finale", EnsembleRank, MyRank, HLTR1HL

            CovSeik1=CovSeik1+HLTR1HL
            TempMatrixSeik=CovSeik1
            call dposv( 'U', SeikDim, 1, TempMatrixSeik, SeikDim, TempSliceSeik2, SeikDim, ierr)
            if (ierr.ne.0) error stop 'Analysis Matrix inversion failed!'            
        else
            call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
        end if
    else
        call mpi_allgatherv(ObsBaseMember,ObsSpaceDim,mpi_real8,HLSeik,MpiCountObs,MpiDisplacementObs,mpi_real8,EnsembleComm,ierr)
        TempObsSeik=reshape(ObsBaseMember,(/ObsSpaceDim/))
        TempObsSeik=TempObsSeik*ObsErrorDiag1
        TempSliceSeik=matmul(TempObsSeik,HLSeik)

if (TempSliceSeik(1)>1) write(*,*) "controllo", EnsembleRank, MyRank, TempSliceSeik
if (EnsembleRank==1) then
    open(unit=UnitSEIK, file="temp"//int2str(MyRank,3)//".csv", form='formatted', iostat=ierr, action='write', access='sequential',status='replace')
    write(UnitSEIK,*,iostat=ierr) ObsBaseMember
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    do jn=1,SeikDim
        write(UnitSEIK,*,iostat=ierr) "----------------------",jn,"---------------------------"
        write(UnitSEIK,*,iostat=ierr) HLSeik(:,jn)
    end do
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) TempSliceSeik
    close(unit=UnitSEIK, iostat=ierr)
end if

        if (MyRank==0) then
            call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,0,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
        else
            call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
        end if
    end if
    
    if (MyRank==0) then
        ChangeCoefSeik=TempSliceSeik2
        call MPI_Bcast( ChangeCoefSeik, SeikDim, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
    end if
    
    call MPI_Bcast( ChangeCoefSeik, SeikDim, mpi_real8, 0, LocalComm, ierr)
    TempVecSeik=matmul(Lseik,ChangeCoefSeik)
    BaseMember=reshape(TempVecSeik,(/ jpk,jpj,jpi,jptra /))
    trnEnsembleWeighted=BaseMember**2
    where (trnEnsembleWeighted>MaxVarSEIK) BaseMember=sign(sqrt(MaxVarSEIK),BaseMember)
    
    BaseMember=exp(BaseMember)
    trn=trn*BaseMember
    !trb=trn
end subroutine
