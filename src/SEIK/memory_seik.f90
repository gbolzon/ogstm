      MODULE memory_seik
       
      USE modul_param

      IMPLICIT NONE

      public

      INTEGER :: SeikDim ! Reduced dimension of the error subspace.
      integer :: SpaceDim ! Full dimension of the space
      INTEGER :: EnsembleComm, EnsembleRank, EnsembleSize !, BaseComm
      
      integer, parameter :: NotWorkingMember=0, UnitSEIK=1001
      logical, parameter :: UseInflation=.false., UseHighOrder=.false., UseModSeik=.false., UseMaxVarSEIK=.true., UseDiffCov=.true.
      character(len=*), parameter :: PCANeeded="none" ! "read" = read the matrices in the SAVE folder and do pca, "write"= save the matrices and do pca, anything else means no pca 
      logical, parameter :: PCAFullYear=.false.
      double precision, parameter :: MaxVarSEIK=1.0d0, CutOffValue=1.0d-5

      double precision, allocatable, dimension (:,:,:,:) :: trnEnsemble, trnEnsembleWeighted, BaseMember
      double precision, allocatable, dimension (:) :: ModelErrorDiag1
      double precision, allocatable, dimension (:,:) :: LSeik
      double precision :: Weight ! teporary defined in myalloc_seik
      double precision, allocatable, dimension (:,:) :: ChangeBaseSeik
      double precision, allocatable, dimension (:) :: ChangeCoefSeik
      
      integer, allocatable, dimension (:) :: MpiCount, MpiDisplacement, MpiCountCov, MpiDisplacementCov
      double precision, allocatable, dimension (:) :: TempVecSeik
      double precision, allocatable, dimension (:) :: TempSliceSeik, TempSliceSeik2
      double precision, allocatable, dimension (:,:) :: TempMatrixSeik, OrtMatrixSampling
      
      !for 1 process
      double precision, allocatable, dimension(:,:) :: CovSeik1, TTTSeik, CovSmoother1part, LTQ1L
      double precision, allocatable, dimension (:) :: AllWeights, AllWeightsSqrt, AllWeightsSqrt1
      double precision :: ForgettingFactor
      
      !for PCANeeded
      integer :: DimForPCA, nHistoryForVar, nHistoryForSVD, nHistoryForSVDpart, CounterForVar, CounterForSVD, CounterForSVDpart, SVDpartID
      double precision, allocatable, dimension(:,:,:,:,:) :: HistoryForVar, HistoryForSVD, HistoryForSVDpart
      double precision, allocatable, dimension(:,:) :: PCAMatrix
      double precision, allocatable, dimension(:,:) :: PCAMatrixT

      double precision,allocatable,dimension(:,:,:) :: copy_inSeik

      double precision, allocatable, dimension (:,:,:,:) :: trnVariance
      
      !for analisis
      integer :: ObsSpaceDim
      double precision, allocatable, dimension(:,:) :: ObsDataSeik, ComputedObsSeik, MisfitSeik
      double precision, allocatable, dimension(:) :: ObsErrorDiag1
      double precision, allocatable, dimension(:,:) :: ComputedObsMean, ObsBaseMember
      double precision, allocatable, dimension(:,:) :: HLSeik 
      integer, allocatable, dimension (:) :: MpiCountObs, MpiDisplacementObs
      double precision, allocatable, dimension (:) :: TempObsSeik
      double precision, allocatable, dimension(:,:) :: HLTR1HL
      
      logical :: CutLeft, CutRight, CutTop, CutBottom ! definied in myalloc
      
      double precision,  dimension(:,:), allocatable :: BufferMPILinkSend3, BufferMPILinkSend4, BufferMPILinkReceive3, BufferMPILinkReceive4
      
      integer, dimension(:,:,:), allocatable :: SeikMask
      
      ! Per UseDiffCov
      integer :: UDiffSpaceDim, VDiffSpaceDim, WDiffSpaceDim, UDiffObsSpaceDim, VDiffObsSpaceDim
      integer, dimension(:,:,:), allocatable :: SeikUMask, SeikVMask, SeikWMask
      double precision, allocatable, dimension (:) :: UDiffModelErrorDiag1, VDiffModelErrorDiag1, WDiffModelErrorDiag1, UDiffObsErrorDiag1, VDiffObsErrorDiag1
      double precision, allocatable, dimension (:,:,:,:) :: UDiffBaseMember, VDiffBaseMember, WDiffBaseMember
      double precision, allocatable, dimension (:,:,:,:) :: TempUDiffBaseMember, TempVDiffBaseMember, TempWDiffBaseMember
      double precision, allocatable, dimension (:) :: TempUDiffVecSeik, TempVDiffVecSeik, TempWDiffVecSeik
      double precision, allocatable, dimension(:,:) :: UDiffObsBaseMember, VDiffObsBaseMember
      double precision, allocatable, dimension(:) :: TempUDiffObsSeik, TempVDiffObsSeik
      double precision, allocatable, dimension (:,:) :: ULSeik, VLSeik, WLSeik, UHLSeik, VHLSeik
      double precision, allocatable, dimension (:,:,:,:) :: UVariance, VVariance, WVariance
      integer, allocatable, dimension (:) :: UDiffMpiCount, UDiffMpiDisplacement,VDiffMpiCount, VDiffMpiDisplacement,WDiffMpiCount, WDiffMpiDisplacement
      integer, allocatable, dimension (:) :: UDiffMpiCountObs, UDiffMpiDisplacementObs, VDiffMpiCountObs, VDiffMpiDisplacementObs
      
      
      CONTAINS
       
      subroutine myalloc_seik(LocalRank)
            implicit none

            integer, intent(in) :: LocalRank        
            integer :: indexi

            Weight=1.0d0/(SeikDim+1) ! it needs a better initialization after AllWeights
            
            SpaceDim=jpk*jpj*jpi*jptra
            
            ObsSpaceDim=jpj*jpi
            
            ForgettingFactor=1.0d0 !0.9d0
            
            !CutLeft=(.not.(nldi==1)) !these are definied in myalloc, here we don't have nldi etc...
            !CutRight=(.not.(nlei==jpi))
            !CutTop=(.not.(nlej==jpj))
            !CutBottom=(.not.(nldj==1))
            
            allocate(ModelErrorDiag1(SpaceDim))
            ModelErrorDiag1 = huge(ModelErrorDiag1(1))
            ModelErrorDiag1 = 1/(log(1.1d0)**2)*1000 !il fattore mille significa che stiamo considerando la varianza sulla superficie
            
            allocate(ObsErrorDiag1(ObsSpaceDim))                    
            ObsErrorDiag1 = huge(ObsErrorDiag1(1))
            ObsErrorDiag1 = 1/(log(1.1d0)**2)
            
            if (UseDiffCov) then
                UDiffSpaceDim=jpk*jpj*(jpi-1)*jptra
                VDiffSpaceDim=jpk*(jpj-1)*jpi*jptra
                WDiffSpaceDim=(jpk-1)*jpj*jpi*jptra
                UDiffObsSpaceDim=jpj*(jpi-1)
                VDiffObsSpaceDim=(jpj-1)*jpi
                
                allocate(UDiffModelErrorDiag1(UDiffSpaceDim))
                UDiffModelErrorDiag1 = huge(UDiffModelErrorDiag1(1))
                UDiffModelErrorDiag1 = 1/(log(1.1d0)**2)*3000 !3000 e' un'approssimazione della massima profondita'. e' 3 volte mille, perche' la prima cella ha profondita' 3
                
                allocate(VDiffModelErrorDiag1(VDiffSpaceDim))
                VDiffModelErrorDiag1 = huge(VDiffModelErrorDiag1(1))
                VDiffModelErrorDiag1 = 1/(log(1.1d0)**2)*3000
                
                allocate(WDiffModelErrorDiag1(WDiffSpaceDim))
                WDiffModelErrorDiag1 = huge(WDiffModelErrorDiag1(1))
                WDiffModelErrorDiag1 = 0.0d0
                
                allocate(UDiffObsErrorDiag1(UDiffObsSpaceDim))
                UDiffObsErrorDiag1 = huge(UDiffObsErrorDiag1(1))
                UDiffObsErrorDiag1 = 1/(log(1.1d0)**2)
                
                allocate(VDiffObsErrorDiag1(VDiffObsSpaceDim))
                VDiffObsErrorDiag1 = huge(VDiffObsErrorDiag1(1))
                VDiffObsErrorDiag1 = 1/(log(1.1d0)**2)
                
                allocate(SeikUMask(jpk,jpj,jpi-1))
                SeikUMask = huge(SeikUMask(1,1,1))
                
                allocate(SeikVMask(jpk,jpj-1,jpi))
                SeikVMask = huge(SeikVMask(1,1,1))
                
                allocate(SeikWMask(jpk-1,jpj,jpi))
                SeikWMask = huge(SeikWMask(1,1,1))
                
                allocate(UDiffBaseMember(jpk,jpj,jpi-1,jptra))
                UDiffBaseMember = huge(UDiffBaseMember(1,1,1,1))
                
                allocate(VDiffBaseMember(jpk,jpj-1,jpi,jptra))
                VDiffBaseMember = huge(VDiffBaseMember(1,1,1,1))
                
                allocate(WDiffBaseMember(jpk-1,jpj,jpi,jptra))
                WDiffBaseMember = huge(WDiffBaseMember(1,1,1,1))
                
                allocate(TempUDiffBaseMember(jpk,jpj,jpi-1,jptra))
                TempUDiffBaseMember = huge(TempUDiffBaseMember(1,1,1,1))
                
                allocate(TempVDiffBaseMember(jpk,jpj-1,jpi,jptra))
                TempVDiffBaseMember = huge(TempVDiffBaseMember(1,1,1,1))
                
                allocate(TempWDiffBaseMember(jpk-1,jpj,jpi,jptra))
                TempWDiffBaseMember = huge(TempWDiffBaseMember(1,1,1,1))
                
                allocate(UVariance(jpk,jpj,jpi-1,jptra))
                UVariance = huge(UVariance(1,1,1,1))
                
                allocate(VVariance(jpk,jpj-1,jpi,jptra))
                VVariance = huge(VVariance(1,1,1,1))
                
                allocate(WVariance(jpk-1,jpj,jpi,jptra))
                WVariance = huge(WVariance(1,1,1,1))
                
                allocate(UDiffObsBaseMember(jpj,jpi-1))                    
                UDiffObsBaseMember = huge(UDiffObsBaseMember(1,1))
                
                allocate(VDiffObsBaseMember(jpj-1,jpi))                    
                VDiffObsBaseMember = huge(VDiffObsBaseMember(1,1))
                
                allocate(TempUDiffObsSeik(UDiffObsSpaceDim))                    
                TempUDiffObsSeik = huge(TempUDiffObsSeik(1))
                
                allocate(TempVDiffObsSeik(VDiffObsSpaceDim))                    
                TempVDiffObsSeik = huge(TempVDiffObsSeik(1))
                
                allocate(ULSeik(UDiffSpaceDim,SeikDim))
                ULSeik = huge(ULSeik(1,1))
                
                allocate(VLSeik(VDiffSpaceDim,SeikDim))
                VLSeik = huge(VLSeik(1,1))
                
                allocate(WLSeik(WDiffSpaceDim,SeikDim))
                WLSeik = huge(WLSeik(1,1))
                
                allocate(UHLSeik(UDiffObsSpaceDim,SeikDim))                    
                UHLSeik = huge(UHLSeik(1,1))
                
                allocate(VHLSeik(VDiffObsSpaceDim,SeikDim))                    
                VHLSeik = huge(VHLSeik(1,1))        
                
                allocate(TempUDiffVecSeik(UDiffSpaceDim))
                TempUDiffVecSeik = huge(TempUDiffVecSeik(1))
                
                allocate(TempVDiffVecSeik(VDiffSpaceDim))
                TempVDiffVecSeik = huge(TempVDiffVecSeik(1))
                
                allocate(TempWDiffVecSeik(WDiffSpaceDim))
                TempWDiffVecSeik = huge(TempWDiffVecSeik(1))
                
                allocate(UDiffMpiCount(0:SeikDim))
                UDiffMpiCount = huge(UDiffMpiCount(1))
                UDiffMpiCount=UDiffSpaceDim
                UDiffMpiCount(NotWorkingMember)=0
                
                allocate(UDiffMpiDisplacement(0:SeikDim))
                UDiffMpiDisplacement = huge(UDiffMpiDisplacement(1))
                UDiffMpiDisplacement(0)=0
                do indexi=1, SeikDim
                    UDiffMpiDisplacement(indexi)=UDiffMpiDisplacement(indexi-1)+UDiffMpiCount(indexi-1)
                end do
                
                allocate(VDiffMpiCount(0:SeikDim))
                VDiffMpiCount = huge(VDiffMpiCount(1))
                VDiffMpiCount=VDiffSpaceDim
                VDiffMpiCount(NotWorkingMember)=0
                
                allocate(VDiffMpiDisplacement(0:SeikDim))
                VDiffMpiDisplacement = huge(VDiffMpiDisplacement(1))
                VDiffMpiDisplacement(0)=0
                do indexi=1, SeikDim
                    VDiffMpiDisplacement(indexi)=VDiffMpiDisplacement(indexi-1)+VDiffMpiCount(indexi-1)
                end do
                
                allocate(WDiffMpiCount(0:SeikDim))
                WDiffMpiCount = huge(WDiffMpiCount(1))
                WDiffMpiCount=WDiffSpaceDim
                WDiffMpiCount(NotWorkingMember)=0
                
                allocate(WDiffMpiDisplacement(0:SeikDim))
                WDiffMpiDisplacement = huge(WDiffMpiDisplacement(1))
                WDiffMpiDisplacement(0)=0
                do indexi=1, SeikDim
                    WDiffMpiDisplacement(indexi)=WDiffMpiDisplacement(indexi-1)+WDiffMpiCount(indexi-1)
                end do
                
                allocate(UDiffMpiCountObs(0:SeikDim))
                UDiffMpiCountObs = huge(UDiffMpiCountObs(1))
                UDiffMpiCountObs=UDiffObsSpaceDim
                UDiffMpiCountObs(NotWorkingMember)=0
                
                allocate(UDiffMpiDisplacementObs(0:SeikDim))
                UDiffMpiDisplacementObs = huge(UDiffMpiDisplacementObs(1))
                UDiffMpiDisplacementObs(0)=0
                do indexi=1, SeikDim
                    UDiffMpiDisplacementObs(indexi)=UDiffMpiDisplacementObs(indexi-1)+UDiffMpiCountObs(indexi-1)
                end do
                
                allocate(VDiffMpiCountObs(0:SeikDim))
                VDiffMpiCountObs = huge(VDiffMpiCountObs(1))
                VDiffMpiCountObs=VDiffObsSpaceDim
                VDiffMpiCountObs(NotWorkingMember)=0
                
                allocate(VDiffMpiDisplacementObs(0:SeikDim))
                VDiffMpiDisplacementObs = huge(VDiffMpiDisplacementObs(1))
                VDiffMpiDisplacementObs(0)=0
                do indexi=1, SeikDim
                    VDiffMpiDisplacementObs(indexi)=VDiffMpiDisplacementObs(indexi-1)+VDiffMpiCountObs(indexi-1)
                end do
                
            end if
            
            allocate(trnEnsemble(jpk,jpj,jpi,jptra))                    
            trnEnsemble = huge(trnEnsemble(1,1,1,1))
            
            allocate(trnEnsembleWeighted(jpk,jpj,jpi,jptra))
            trnEnsembleWeighted = huge(trnEnsembleWeighted(1,1,1,1))
            
            allocate(BaseMember(jpk,jpj,jpi,jptra))
            BaseMember = huge(BaseMember(1,1,1,1))
            
            allocate(LSeik(SpaceDim,SeikDim))
            LSeik = huge(LSeik(1,1))
            
            allocate(MpiCount(0:SeikDim))
            MpiCount = huge(MpiCount(1))
            MpiCount=SpaceDim
            MpiCount(NotWorkingMember)=0
            
            allocate(MpiDisplacement(0:SeikDim))
            MpiDisplacement = huge(MpiDisplacement(1))
            MpiDisplacement(0)=0
            do indexi=1, SeikDim
                MpiDisplacement(indexi)=MpiDisplacement(indexi-1)+MpiCount(indexi-1)
            end do
            
            allocate(TempVecSeik(SpaceDim))
            TempVecSeik = huge(TempVecSeik(1))
            
            allocate(TempSliceSeik(SeikDim))
            TempSliceSeik = huge(TempSliceSeik(1))
            
            allocate(TempMatrixSeik(SeikDim,SeikDim))
            TempMatrixSeik = huge(TempMatrixSeik(1,1))
            
            allocate(ChangeCoefSeik(SeikDim))
            ChangeCoefSeik = huge(ChangeCoefSeik(1))
            
            allocate(ObsDataSeik(jpj,jpi))                    
            ObsDataSeik = huge(ObsDataSeik(1,1))
            
            allocate(ComputedObsSeik(jpj,jpi))                    
            ComputedObsSeik = huge(ComputedObsSeik(1,1))
            
            allocate(MisfitSeik(jpj,jpi))                    
            MisfitSeik = huge(MisfitSeik(1,1))
            
            allocate(ComputedObsMean(jpj,jpi))                    
            ComputedObsMean = huge(ComputedObsMean(1,1))
            
            allocate(ObsBaseMember(jpj,jpi))                    
            ObsBaseMember = huge(ObsBaseMember(1,1))
            
            allocate(HLSeik(ObsSpaceDim,SeikDim))                    
            HLSeik = huge(HLSeik(1,1))
            
            allocate(MpiCountObs(0:SeikDim))
            MpiCountObs = huge(MpiCountObs(1))
            MpiCountObs=ObsSpaceDim
            MpiCountObs(NotWorkingMember)=0
            
            allocate(MpiDisplacementObs(0:SeikDim))
            MpiDisplacementObs = huge(MpiDisplacementObs(1))
            MpiDisplacementObs(0)=0
            do indexi=1, SeikDim
                MpiDisplacementObs(indexi)=MpiDisplacementObs(indexi-1)+MpiCountObs(indexi-1)
            end do
            
            allocate(TempObsSeik(ObsSpaceDim))                    
            TempObsSeik = huge(TempObsSeik(1))
            
            allocate(BufferMPILinkSend3(jpk,jpi))                    
            BufferMPILinkSend3 = huge(BufferMPILinkSend3(1,1))
            
            allocate(BufferMPILinkSend4(jpk,jpi))                    
            BufferMPILinkSend4 = huge(BufferMPILinkSend4(1,1))
            
            allocate(BufferMPILinkReceive3(jpk,jpi))                    
            BufferMPILinkReceive3 = huge(BufferMPILinkReceive3(1,1))
            
            allocate(BufferMPILinkReceive4(jpk,jpi))                    
            BufferMPILinkReceive4 = huge(BufferMPILinkReceive4(1,1))
            
            allocate(SeikMask(jpk,jpj,jpi))
            SeikMask = huge(SeikMask(1,1,1))

            allocate(trnVariance(jpk,jpj,jpi,jptra))
            trnVariance = huge(trnVariance(1,1,1,1))

            if(LocalRank==0) then
                
                allocate(TempSliceSeik2(SeikDim))
                TempSliceSeik2 = huge(TempSliceSeik2(1))
                
                allocate(MpiCountCov(0:SeikDim))
                MpiCountCov = huge(MpiCountCov(1))
                MpiCountCov=SeikDim
                MpiCountCov(NotWorkingMember)=0
                
                allocate(MpiDisplacementCov(0:SeikDim))
                MpiDisplacementCov = huge(MpiDisplacementCov(1))
                MpiDisplacementCov(0)=0
                do indexi=1, SeikDim
                    MpiDisplacementCov(indexi)=MpiDisplacementCov(indexi-1)+MpiCountCov(indexi-1)
                end do
                
                allocate(copy_inSeik(jpiglo, jpjglo, jpk))
                copy_inSeik = huge(copy_inSeik(1,1,1))
                    
            end if
            
            
            if(EnsembleRank==NotWorkingMember) then
                
                allocate(AllWeights(0:SeikDim))                    
                AllWeights = huge(AllWeights(0))
                AllWeights = 1.0d0/(SeikDim+1)
                
                allocate(TTTSeik(SeikDim,SeikDim))
                TTTSeik = huge(TTTSeik(1,1))
                call TTTSeik_builder()
                
                
                if(LocalRank==0) then
                
                    allocate(AllWeightsSqrt(0:SeikDim))                    
                    AllWeightsSqrt = huge(AllWeightsSqrt(0))
                    AllWeightsSqrt = sqrt(AllWeights)
                
                    allocate(AllWeightsSqrt1(0:SeikDim))                    
                    AllWeightsSqrt1 = huge(AllWeightsSqrt1(0))
                    AllWeightsSqrt1 = 1.0d0/AllWeightsSqrt
            
                    allocate(CovSeik1(SeikDim,SeikDim))
                    CovSeik1 = huge(CovSeik1(1,1))
                    
                    allocate(CovSmoother1part(SeikDim,SeikDim))
                    CovSmoother1part = huge(CovSmoother1part(1,1))
                    
                    allocate(LTQ1L(SeikDim,SeikDim))
                    LTQ1L = huge(LTQ1L(1,1))
                    
                    allocate(OrtMatrixSampling(SeikDim+1,SeikDim+1))
                    OrtMatrixSampling = huge(OrtMatrixSampling(1,1))
                    
                    allocate(ChangeBaseSeik(SeikDim,0:SeikDim))
                    ChangeBaseSeik = huge(ChangeBaseSeik(1,1))
                    
                    
                    allocate(HLTR1HL(SeikDim,SeikDim))
                    HLTR1HL = huge(HLTR1HL(1,1))
                    
                end if
                
            end if
            
            if ((PCANeeded.eq."read").or.(PCANeeded.eq."write")) then
            
                DimForPCA=100
                nHistoryForVar=365
                nHistoryForSVD=31*10
                nHistoryForSVDpart=31
                CounterForVar=0
                CounterForSVD=0 
                CounterForSVDpart=0 
                SVDpartID=0
                
                if (PCAFullYear) then
                    allocate(HistoryForVar(jpk,jpj,jpi,jptra,nHistoryForVar))
                    HistoryForVar = huge(HistoryForVar(1,1,1,1,1))
                end if
        
                allocate(HistoryForSVD(jpk,jpj,jpi,jptra,nHistoryForSVD))
                HistoryForSVD = huge(HistoryForSVD(1,1,1,1,1))
                
                allocate(HistoryForSVDpart(jpk,jpj,jpi,jptra,nHistoryForSVDpart))
                HistoryForSVDpart = huge(HistoryForSVDpart(1,1,1,1,1))

                allocate(PCAMatrix(nHistoryForSVD,SpaceDim))
                PCAMatrix = huge(PCAMatrix(1,1))

                allocate(PCAMatrixT(SpaceDim,nHistoryForSVD))
                PCAMatrixT = huge(PCAMatrixT(1,1))

            end if
            
      end subroutine 
      
      subroutine TTTSeik_builder()
            implicit none
            double precision, dimension(0:SeikDim,SeikDim) :: matrixT
            integer :: indexi
            
            do indexi=1, SeikDim
                matrixT(:,indexi)=-AllWeights
                matrixT(indexi-1,indexi)=matrixT(indexi-1,indexi)+1
            end do
            do indexi=1, SeikDim
                TTTSeik(:,indexi)=matmul(matrixT(:,indexi)/AllWeights,matrixT)
            end do
            
      end subroutine

      subroutine clean_seik(LocalRank)
            implicit none

            integer, intent(in) :: LocalRank

            deallocate(trnEnsemble)
            deallocate(trnEnsembleWeighted)
            deallocate(BaseMember)
            deallocate(ModelErrorDiag1)
            deallocate(LSeik)
            deallocate(MpiCount)
            deallocate(MpiDisplacement)
            deallocate(TempVecSeik)
            deallocate(TempSliceSeik)
            deallocate(TempMatrixSeik)
            deallocate(ChangeCoefSeik)
            deallocate(ObsDataSeik)
            deallocate(ComputedObsSeik)
            deallocate(MisfitSeik)
            deallocate(ObsErrorDiag1)
            deallocate(ComputedObsMean)
            deallocate(ObsBaseMember)
            deallocate(HLSeik)
            deallocate(MpiCountObs)
            deallocate(MpiDisplacementObs)
            deallocate(TempObsSeik)
            deallocate(BufferMPILinkSend3)                    
            deallocate(BufferMPILinkSend4)                    
            deallocate(BufferMPILinkReceive3)                    
            deallocate(BufferMPILinkReceive4)                    
            deallocate(SeikMask)
            deallocate(trnVariance)
            
            if (UseDiffCov) then                
                deallocate(UDiffModelErrorDiag1)
                deallocate(VDiffModelErrorDiag1)
                deallocate(WDiffModelErrorDiag1)
                deallocate(UDiffObsErrorDiag1)
                deallocate(VDiffObsErrorDiag1)
                deallocate(SeikUMask)
                deallocate(SeikVMask)
                deallocate(SeikWMask)
                deallocate(UDiffBaseMember)
                deallocate(VDiffBaseMember)
                deallocate(WDiffBaseMember)
                deallocate(TempUDiffBaseMember)
                deallocate(TempVDiffBaseMember)
                deallocate(TempWDiffBaseMember)
                deallocate(UVariance)
                deallocate(VVariance)
                deallocate(WVariance)
                deallocate(UDiffObsBaseMember)                    
                deallocate(VDiffObsBaseMember)                    
                deallocate(ULSeik)
                deallocate(VLSeik)
                deallocate(WLSeik)
                deallocate(UHLSeik)                    
                deallocate(VHLSeik)                    
                deallocate(TempUDiffVecSeik)
                deallocate(TempVDiffVecSeik)
                deallocate(TempWDiffVecSeik)
                deallocate(TempUDiffObsSeik)
                deallocate(TempVDiffObsSeik)
                deallocate(UDiffMpiCount)
                deallocate(UDiffMpiDisplacement)
                deallocate(VDiffMpiCount)
                deallocate(VDiffMpiDisplacement)
                deallocate(WDiffMpiCount)
                deallocate(WDiffMpiDisplacement)
                deallocate(UDiffMpiCountObs)
                deallocate(UDiffMpiDisplacementObs)
                deallocate(VDiffMpiCountObs)
                deallocate(VDiffMpiDisplacementObs)
                
            end if
           
            if (LocalRank==0) then
                deallocate(TempSliceSeik2)
                deallocate(MpiCountCov)
                deallocate(MpiDisplacementCov)
                deallocate(copy_inSeik)
            end if
            
            if (EnsembleRank==0) then
                deallocate(AllWeights)
                deallocate(TTTSeik)
                
                if (LocalRank==0) then
                    deallocate(AllWeightsSqrt)
                    deallocate(AllWeightsSqrt1)
                    deallocate(CovSeik1)
                    deallocate(CovSmoother1part)
                    deallocate(LTQ1L)
                    deallocate(OrtMatrixSampling)
                    deallocate(ChangeBaseSeik)
                    deallocate(HLTR1HL)
                end if
            end if
            
            if ((PCANeeded.eq."read").or.(PCANeeded.eq."write")) then
                if (PCAFullYear) then                
                    deallocate(HistoryForVar)
                end if                
                deallocate(HistoryForSVD)
                deallocate(HistoryForSVDpart)
                deallocate(PCAMatrix)
                deallocate(PCAMatrixT)
            end if
       
      end subroutine

      END MODULE memory_seik
