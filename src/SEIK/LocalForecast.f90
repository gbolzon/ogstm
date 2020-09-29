subroutine LocalForecast
    Use myalloc
    use mpi


    implicit none
    INTEGER :: ierr, indexj, indexi, indexk, indexl, temp

    trnEnsemble=trn
    where (trnEnsemble<CutOffValue) 
        BaseMember=log(CutOffValue)
    elsewhere
        BaseMember=log(trnEnsemble)
    end where
    trnEnsembleWeighted=BaseMember*SeikWeight
    
    call MPI_AllReduce(trnEnsembleWeighted, trn, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    
    BaseMember=BaseMember-trn
    where (SeikTrcMask==0) BaseMember=0.0d0   
    
    do indexi=1, jptra
        where (tmask==1) 
            trnEnsembleWeighted(:,:,:,indexi)=log(trnEnsemble(:,:,:,indexi))*SeikWeight
        elsewhere
            trnEnsembleWeighted(:,:,:,indexi)=0.0d0
        end where
    end do
    
    call MPI_AllReduce(trnEnsembleWeighted, trn, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)

    !if (CutLeft) BaseMember(:,:,1,:)=0.0d0
    !if (CutRight) BaseMember(:,:,jpi,:)=0.0d0
    !if (CutTop) BaseMember(:,jpj,:,:)=0.0d0
    !if (CutBottom) BaseMember(:,1,:,:)=0.0d0

    if (UseDiffCov) then !to be done ,not working atm
        do indexi=1, jptra
            do indexj=1, jpk
                UDiffBaseMember(indexj,:,:,indexi)=(BaseMember(indexj,:,2:jpi,indexi)-BaseMember(indexj,:,1:jpi-1,indexi))/e1u(:,1:jpi-1)
            end do
            UDiffBaseMember(:,:,:,indexi)=UDiffBaseMember(:,:,:,indexi)*SeikUMask
        end do
        
        do indexi=1, jptra
            do indexj=1, jpk
                VDiffBaseMember(indexj,:,:,indexi)=(BaseMember(indexj,2:jpj,:,indexi)-BaseMember(indexj,1:jpj-1,:,indexi))/e2v(1:jpj-1,:)
            end do
            VDiffBaseMember(:,:,:,indexi)=VDiffBaseMember(:,:,:,indexi)*SeikVMask
        end do
        
        do indexi=1, jptra
            WDiffBaseMember(:,:,:,indexi)=(BaseMember(2:jpk,:,:,indexi)-BaseMember(1:jpk-1,:,:,indexi))/e3w_0(2:jpk,:,:)
            WDiffBaseMember(:,:,:,indexi)=WDiffBaseMember(:,:,:,indexi)*SeikWMask
        end do
        
    end if
    
    trnEnsembleWeighted=BaseMember**2*SeikWeight
    call MPI_AllReduce(trnEnsembleWeighted, trnVariance, SpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
    
    if (UseMaxVarSEIK) then
        where (trnVariance>MaxVarVec) BaseMember=BaseMember*sqrt(MaxVarVec/trnVariance)
        
        !qui bisognerebbe introdurre maxvarUvec ecc
        if (UseDiffCov) then !to be done, not working atm
            TempUDiffBaseMember=UDiffBaseMember**2*SeikWeight
            call MPI_AllReduce(TempUDiffBaseMember, UVariance, UDiffSpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
            where (UVariance>MaxVarSEIK) UDiffBaseMember=UDiffBaseMember*sqrt(MaxVarSEIK/UVariance)
            
            TempVDiffBaseMember=VDiffBaseMember**2*SeikWeight
            call MPI_AllReduce(TempVDiffBaseMember, VVariance, VDiffSpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
            where (VVariance>MaxVarSEIK) VDiffBaseMember=VDiffBaseMember*sqrt(MaxVarSEIK/VVariance)
            
            TempWDiffBaseMember=WDiffBaseMember**2*SeikWeight
            call MPI_AllReduce(TempWDiffBaseMember, WVariance, WDiffSpaceDim, mpi_real8, MPI_SUM, EnsembleComm,ierr)
            where (WVariance>MaxVarSEIK) WDiffBaseMember=WDiffBaseMember*sqrt(MaxVarSEIK/WVariance)
            
        end if
    end if

    if (EnsembleRank==NotWorkingMember) then
        call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        
        if (UseDiffCov) then !to be done, not working atm
            call mpi_allgatherv(0,0,mpi_real8,ULSeik,UDiffMpiCount,UDiffMpiDisplacement,mpi_real8,EnsembleComm,ierr)
            call mpi_allgatherv(0,0,mpi_real8,VLSeik,VDiffMpiCount,VDiffMpiDisplacement,mpi_real8,EnsembleComm,ierr)
            call mpi_allgatherv(0,0,mpi_real8,WLSeik,WDiffMpiCount,WDiffMpiDisplacement,mpi_real8,EnsembleComm,ierr)
        end if
        
        if (UseInflation) then
            CovSeik1=ForgettingFactor*TTTSeik
        else
        
            call mpi_gatherv(0,0,mpi_real8,LTQ1L_sjis,LocalMpiCountCov,LocalMpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
                       
!write(*,*) "myrank=", myrank, " EnsembleRank=", EnsembleRank, " riga 99"
call mpi_barrier(LocalComm, ierr)
if (myRank==40) then
    write(*,*) "---------------------------------------------------------"
    write(*,*) "LTQ1L_sjis"
    write(*,*) "rank 40, x=5, y=5"
    do indexi=1, SeikDim
        write(*,*) LTQ1L_sjis(:,5,5, indexi)
    end do
    write(*,*) "---------------------------------------------------------"
end if
call mpi_barrier(LocalComm, ierr)

            do indexi=1, jpi
                do indexj=1, jpj
                    if (SeikMask(1, indexj,indexi)==1) then
                    
                        LTQ1L=LTQ1L_sjis(:,indexj,indexi,:)
                        
                        call ForecastMatrixOp
                        
                        call SymChangeBase(CovSeik1,indexi, indexj)
                        
                        LTQ1L_sjis(:,indexj,indexi,:)=CovSeik1
                        
                    else
                        LTQ1L_sjis(:,indexj,indexi,:)=0.0d0
                    end if
                end do
            end do
            
!write(*,*) "myrank=", myrank, " EnsembleRank=", EnsembleRank, " riga 117"

            call MPI_Scatterv(LTQ1L_sjis, LocalMpiCountCov, LocalMpiDisplacementCov, mpi_real8, 0, 0, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
            
            CovSeik1=0.0d0
            do indexi=1,SeikDim
                CovSeik1(indexi,indexi)=1.0d0
            end do
            
!write(*,*) "myrank=", myrank, " EnsembleRank=", EnsembleRank, " riga 126"
            
            call mpi_allgatherv(0,0,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)

!write(*,*) "myrank=", myrank, " EnsembleRank=", EnsembleRank, " riga 130"

        end if
                        
if (.false.) then !non locale
        if (MyRank==0) then
            if (UseInflation) then
                CovSeik1=ForgettingFactor*TTTSeik
            else
                call mpi_gatherv(0,0,mpi_real8,LTQ1L,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)

write(*,*) "LTQ1L=", LTQ1L

                call ForecastMatrixOp
            end if
        end if
end if    

    else
        call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
        
        if (UseDiffCov) then ! to be done, not working
            call mpi_allgatherv(UDiffBaseMember,UDiffSpaceDim,mpi_real8,ULSeik,UDiffMpiCount,UDiffMpiDisplacement,mpi_real8,EnsembleComm,ierr)
            call mpi_allgatherv(VDiffBaseMember,VDiffSpaceDim,mpi_real8,VLSeik,VDiffMpiCount,VDiffMpiDisplacement,mpi_real8,EnsembleComm,ierr)
            call mpi_allgatherv(WDiffBaseMember,WDiffSpaceDim,mpi_real8,WLSeik,WDiffMpiCount,WDiffMpiDisplacement,mpi_real8,EnsembleComm,ierr)
        end if
        
        if (.not.UseInflation) then
        
            do indexi=1, SeikDim
                TempVecSeik=reshape(BaseMember,(/SpaceDim/))
                TempVecSeik=TempVecSeik*modelErrorDiag1
                TempVecSeik=TempVecSeik*LSeik(:,indexi)
                trnEnsembleWeighted=reshape(TempVecSeik, (/jpk,jpj,jpi,jptra/))
                BaseMember_jic=sum(trnEnsembleWeighted,1)
                BaseMember_ji=sum(BaseMember_jic,3)
                BaseMember_sji(indexi,:,:)=BaseMember_ji
            end do

            call LocalSendRecive(BaseMember_sji)
            
            call SummingLocalPatch(BaseMember_sji)
                      
!write(*,*) "myrank=", myrank, " EnsembleRank=", EnsembleRank, " riga 191"
            
            call mpi_gatherv(BaseMember_sji,SeikDim*jpj*jpi,mpi_real8,0,LocalMpiCountCov,LocalMpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
                            
            call MPI_Scatterv(0, LocalMpiCountCov, LocalMpiDisplacementCov, mpi_real8, BaseMember_sji, SeikDim*jpj*jpi, mpi_real8, NotWorkingMember, EnsembleComm, ierr)
            
!LSeik_reshape=reshape(LSeik, (/jpk, jpj, jpi, jptra, SeikDim/))
            
            call LocalProduct(BaseMember_sji, BaseMember)
            
!write(*,*) "myrank=", myrank, " EnsembleRank=", EnsembleRank, " riga 208"
            
            call mpi_allgatherv(BaseMember,SpaceDim,mpi_real8,LSeik,MpiCount,MpiDisplacement,mpi_real8,EnsembleComm,ierr)
            
if (.false.) then !codice della versione non locale            
            TempVecSeik=reshape(BaseMember,(/SpaceDim/))
            TempVecSeik=TempVecSeik*ModelErrorDiag1
            TempSliceSeik=matmul(TempVecSeik,LSeik)
            
            if (UseDiffCov) then !to be done, not working atm
                TempUDiffVecSeik=reshape(UDiffBaseMember,(/UDiffSpaceDim/))
                TempUDiffVecSeik=TempUDiffVecSeik*UDiffModelErrorDiag1
                TempSliceSeik=TempSliceSeik+matmul(TempUDiffVecSeik,ULSeik)
                
                TempVDiffVecSeik=reshape(VDiffBaseMember,(/VDiffSpaceDim/))
                TempVDiffVecSeik=TempVDiffVecSeik*VDiffModelErrorDiag1
                TempSliceSeik=TempSliceSeik+matmul(TempVDiffVecSeik,VLSeik)
                
                TempWDiffVecSeik=reshape(WDiffBaseMember,(/WDiffSpaceDim/))
                TempWDiffVecSeik=TempWDiffVecSeik*WDiffModelErrorDiag1
                TempSliceSeik=TempSliceSeik+matmul(TempWDiffVecSeik,WLSeik)
            end if

            if (MyRank==0) then
                call mpi_reduce(TempSliceSeik,TempSliceSeik2, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
                call mpi_gatherv(TempSliceSeik2,SeikDim,mpi_real8,0,MpiCountCov,MpiDisplacementCov,mpi_real8,NotWorkingMember, EnsembleComm, ierr)
            else
                call mpi_reduce(TempSliceSeik,0, SeikDim, mpi_real8, mpi_sum, 0, LocalComm, ierr)
            end if
end if
            
        end if
    end if
        
    !end if
    
    trn=exp(trn)
    do indexi=1, jptra
        trn(:,:,:,indexi)=trn(:,:,:,indexi)*tmask
    end do
    
    if (UseCutOffN) then
        do indexi=1,jptra
            if (ctrcnm(indexi).eq."N1p") then
                where (trn(:,:,:,indexi)>CutOffN1p) trn(:,:,:,indexi)=CutOffN1p
            else if (ctrcnm(indexi).eq."N3n") then
                where (trn(:,:,:,indexi)>CutOffN3n) trn(:,:,:,indexi)=CutOffN3n
            end if
        end do
    end if
    
    trb=trn

end subroutine
 
