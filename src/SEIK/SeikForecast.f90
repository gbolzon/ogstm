subroutine SeikForecast(datestring)
    Use myalloc
    use mpi
use StringSEIK

    implicit none
character(len=*), intent(in) :: datestring
    INTEGER :: ierr, indexi
integer :: indexj, indexk, indexn

    where (trn<1.0d-12) trn=1.0d-12 !this should already be done with SMALL in another routine, but I don't have time to check now
    do indexi=1, jptra
        where (tmask==1) 
            trnEnsemble(:,:,:,indexi)=trn(:,:,:,indexi)
            BaseMember(:,:,:,indexi)=log(trn(:,:,:,indexi))
            trnEnsembleWeighted(:,:,:,indexi)=BaseMember(:,:,:,indexi)*Weight
        elsewhere
            trnEnsemble(:,:,:,indexi)=0.0d0
            BaseMember(:,:,:,indexi)=0.0d0
            trnEnsembleWeighted(:,:,:,indexi)=0.0d0
        end where
    end do
    !trnEnsemble=log(trn)
    !trnEnsembleWeighted=trnEnsemble*Weight

    call MPI_AllReduce(trnEnsembleWeighted, trn, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    !BaseMember=trnEnsemble-trn
    where (trn>log(1.0d-5)) 
        BaseMember=BaseMember-trn
    elsewhere
        BaseMember=0.0d0
    end where
    
    !if (CutLeft) BaseMember(:,:,1,:)=0.0d0
    !if (CutRight) BaseMember(:,:,jpi,:)=0.0d0
    !if (CutTop) BaseMember(:,jpj,:,:)=0.0d0
    !if (CutBottom) BaseMember(:,1,:,:)=0.0d0

do indexi=1, jptra
BaseMember(:,:,:,indexi)=BaseMember(:,:,:,indexi)*SeikMask
end do
    
    if (EnsembleRank==NotWorkingMember) then
        call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        if (MyRank==0) then
            if (UseInflation) then
                CovSeik1=ForgettingFactor*TTTSeik
            else
                call mpi_gatherv(0,0,mpi_real8,LTQ1L,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)

write(*,*) "LTQ1L=", LTQ1L

                CovSmoother1part=TTTSeik+LTQ1L
                CovSeik1=TTTSeik
                TempMatrixSeik=CovSmoother1part
                call InvMatMul(TempMatrixSeik,CovSeik1,SeikDim,ierr)
                TempMatrixSeik=matmul(LTQ1L,CovSeik1)
                CovSeik1=TempMatrixSeik
            end if
        end if
    else
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        if (.not.UseInflation) then
            TempVecSeik=reshape(BaseMember,(/SpaceDim/))
            TempVecSeik=TempVecSeik*ModelErrorDiag1
            TempSliceSeik=matmul(TempVecSeik,LSeik)

if (str2int(datestring(1:8))>20130501) then
if (sum(abs(TempSliceSeik))>1) then
    write(*,*) "controllo", EnsembleRank, MyRank, TempSliceSeik
    do indexn=1, jptra
        do indexi=1, jpi
            do indexj=1, jpj
                do indexk=1, jpk
                    if ((abs(BaseMember(indexk, indexj, indexi, indexn))>0.1).and.(abs(ModelErrorDiag1((indexn-1)*jpi*jpj*jpk+(indexi-1)*jpj*jpk+(indexj-1)*jpk+indexk))>1.0d-12)) write(*,*) EnsembleRank, MyRank, trim(ctrcnm(indexn)), "(",indexk, indexj, indexi, indexn,")", BaseMember(indexk, indexj, indexi, indexn), ModelErrorDiag1((indexn-1)*jpi*jpj*jpk+(indexi-1)*jpj*jpk+(indexj-1)*jpk+indexk)  
                end do
            end do
        end do
    end do
if (.false.) then
    !open(unit=UnitSEIK, file="TEMP/temp"//int2str(EnsembleRank,2)//"."//int2str(MyRank,3)//"."//int2str(int(sum(TempSliceSeik)*1000),10)//".csv", form='formatted', iostat=ierr, action='write', access='sequential',status='replace')
    write(UnitSEIK,*,iostat=ierr) BaseMember
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) TempVecSeik
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) TempSliceSeik
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) bfmmask
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) "----------------------------------------------------"
    write(UnitSEIK,*,iostat=ierr) ModelErrorDiag1
    close(unit=UnitSEIK, iostat=ierr)
end if
end if
end if

            if (MyRank==0) then
                call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
                call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,0,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
            else
                call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            end if
        end if
    end if
        
    !end if
    
    trn=exp(trn)
    do indexi=1, jptra
        trn(:,:,:,indexi)=trn(:,:,:,indexi)*tmask
    end do
    trb=trn

end subroutine
