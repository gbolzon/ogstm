      MODULE memory_seik
       
      USE modul_param

      IMPLICIT NONE

      public

      INTEGER :: SeikDim ! Reduced dimension of the error subspace.
      integer :: SpaceDim ! Full dimension of the space
      INTEGER :: EnsembleComm, EnsembleRank, EnsembleSize !, BaseComm
      
      integer, parameter :: NotWorkingMember=0, UnitSEIK=1001
      logical, parameter :: UseInflation=.false.
      character(len=*), parameter :: PCANeeded="none" ! "read" = read the matrices in the SAVE folder and do pca, "write"= save the matrices and do pca, anything else means no pca 
      double precision, parameter :: MaxVarSEIK=4.0d0

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
      double precision :: ForgettingFactor=0.9
      
      !for PCANeeded
      integer :: DimForPCA, nHistoryForVar, nHistoryForSVD, nHistoryForSVDpart, CounterForVar, CounterForSVD, CounterForSVDpart, SVDpartID
      double precision, allocatable, dimension(:,:,:,:,:) :: HistoryForVar, HistoryForSVD, HistoryForSVDpart
      double precision, allocatable, dimension(:,:) :: PCAMatrix
      double precision, allocatable, dimension(:,:) :: PCAMatrixT

      double precision,allocatable,dimension(:,:,:) :: copy_inSeik
      
      !for analisis
      integer :: ObsSpaceDim
      double precision, allocatable, dimension(:,:) :: ObsDataSeik, ComputedObsSeik, MisfitSeik
      double precision, allocatable, dimension(:) :: ObsErrorDiag1
      double precision, allocatable, dimension(:,:) :: ComputedObsMean, ObsBaseMember
      double precision, allocatable, dimension(:,:) :: HLSeik 
      integer, allocatable, dimension (:) :: MpiCountObs, MpiDisplacementObs
      double precision, allocatable, dimension (:) :: TempObsSeik
      double precision, allocatable, dimension(:,:) :: HLTR1HL
      
      
      CONTAINS
       
      subroutine myalloc_seik(LocalRank)
            implicit none

            integer, intent(in) :: LocalRank        
            integer :: indexi

            Weight=1.0d0/(SeikDim+1) ! it needs a better initialization after AllWeights
            SpaceDim=jpk*jpj*jpi*jptra
            ObsSpaceDim=jpj*jpi
            
            allocate(trnEnsemble(jpk,jpj,jpi,jptra))                    
            trnEnsemble = huge(trnEnsemble(1,1,1,1))
            
            allocate(trnEnsembleWeighted(jpk,jpj,jpi,jptra))
            trnEnsembleWeighted = huge(trnEnsembleWeighted(1,1,1,1))
            
            allocate(BaseMember(jpk,jpj,jpi,jptra))
            BaseMember = huge(BaseMember(1,1,1,1))
            
            allocate(ModelErrorDiag1(SpaceDim))
            ModelErrorDiag1 = huge(ModelErrorDiag1(1))
            ModelErrorDiag1 = 1/(log(2.0d0)**2)
            
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
            
            allocate(ObsErrorDiag1(ObsSpaceDim))                    
            ObsErrorDiag1 = huge(ObsErrorDiag1(1))
            
            allocate(ComputedObsMean(jpj,jpi))                    
            ComputedObsMean = huge(ComputedObsMean(1,1))
            
            allocate(ObsBaseMember(jpj,jpi))                    
            ObsBaseMember = huge(ObsBaseMember(1,1))
            
            allocate(HLSeik(ObsSpaceDim,SeikDim))                    
            HLSeik = huge(HLSeik(1,1))
            
            allocate(MpiCountObs(0:ObsSpaceDim))
            MpiCountObs = huge(MpiCountObs(1))
            MpiCountObs=ObsSpaceDim
            MpiCountObs(NotWorkingMember)=0
            
            allocate(MpiDisplacementObs(0:ObsSpaceDim))
            MpiDisplacementObs = huge(MpiDisplacementObs(1))
            MpiDisplacementObs(0)=0
            do indexi=1, ObsSpaceDim
                MpiDisplacementObs(indexi)=MpiDisplacementObs(indexi-1)+MpiCountObs(indexi-1)
            end do
            
            allocate(TempObsSeik(ObsSpaceDim))                    
            TempObsSeik = huge(TempObsSeik(1))
            
            
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
                nHistoryForSVD=310
                nHistoryForSVDpart=31
                CounterForVar=0
                CounterForSVD=0 
                CounterForSVDpart=0 
                SVDpartID=0
            
                !allocate(HistoryForVar(jpk,jpj,jpi,jptra,nHistoryForVar))
                !HistoryForVar = huge(HistoryForVar(1,1,1,1,1))
                
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
                !deallocate(HistoryForVar)
                deallocate(HistoryForSVD)
                deallocate(HistoryForSVDpart)
                deallocate(PCAMatrix)
                deallocate(PCAMatrixT)
            end if
       
      end subroutine

      END MODULE memory_seik
