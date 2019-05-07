      MODULE memory_seik
       
      USE modul_param

      IMPLICIT NONE

      public

      INTEGER :: SeikDim ! Reduced dimension of the error subspace.
      integer :: SpaceDim ! Full dimension of the space
      INTEGER :: EnsembleComm, EnsembleRank, EnsembleSize !, BaseComm
      integer :: NotWorkingMember=0
      logical :: UseInflation=.true.
      double precision, allocatable, dimension (:,:,:,:) :: trnEnsemble, trnEnsembleWeighted, BaseMember
      double precision, allocatable, dimension (:) :: ModelErrorDiag1
      double precision, allocatable, dimension (:,:) :: LSeik
      double precision :: Weight ! teporary defined in myalloc_seik
      
      integer, allocatable, dimension (:) :: MpiCount, MpiDisplacement, MpiCountSeik, MpiDisplacementSeik
      double precision, allocatable, dimension (:) :: TempVecSeik
      double precision, allocatable, dimension (:) :: TempSliceSeik, TempSliceSeik2
      double precision, allocatable, dimension (:,:) :: TempMatrixSeik
      
      !for 1 process
      double precision, allocatable, dimension(:,:) :: CovSeik1, TTTSeik, CovSmoother1part
      double precision, allocatable, dimension (:) :: AllWeights
      double precision :: ForgettingFactor=0.9
       
      CONTAINS
       
      subroutine myalloc_seik()
            Weight=1.0d0/(SeikDim+1) ! it needs a better initialization after AllWeights
            SpaceDim=jpk*jpj*jpi*jptra
            
            allocate(trnEnsemble(jpk,jpj,jpi,jptra))                    
            trnEnsemble = huge(trnEnsemble(1,1,1,1))
            
            allocate(trnEnsembleWeighted(jpk,jpj,jpi,jptra))
            trnEnsembleWeighted = huge(trnEnsembleWeighted(1,1,1,1))
            
            allocate(BaseMember(jpk,jpj,jpi,jptra))
            BaseMember = huge(BaseMember(1,1,1,1))
            
            allocate(ModelErrorDiag1(SpaceDim))
            ModelErrorDiag1 = huge(ModelErrorDiag1(1))
            ModelErrorDiag1 = 1/(log(1.1d0)**2)
            
            allocate(LSeik(SpaceDim,SeikDim))
            LSeik = huge(LSeik(1,1))
            
            allocate(MpiCount(0:SeikDim))
            MpiCount = huge(MpiCount(1))
            MpiCount=SpaceDim
            MpiCount(NotWorkingMember)=0
            
            allocate(MpiDisplacement(0:SeikDim))
            MpiDisplacement = huge(MpiDisplacement(1))
            MpiDisplacement(0)=0
            do i=1, SeikDim
                MpiDisplacement(i)=MpiDisplacement(i-1)+MpiCount(i-1)
            end do
            
            
            
            allocate(TempVecSeik(SpaceDim))
            TempVecSeik = huge(TempVecSeik(1))
            
            allocate(TempSliceSeik(SeikDim))
            TempSliceSeik = huge(TempSliceSeik(1))
            
            allocate(TempMatrixSeik(SeikDim,SeikDim))
            TempMatrixSeik = huge(TempMatrixSeik(1,1))
            
            
            if(MyRank==0) then
                
                allocate(TempSliceSeik2(SeikDim))
                TempSliceSeik2 = huge(TempSliceSeik2(1))
                
                allocate(MpiCountSeik(0:SeikDim))
                MpiCountSeik = huge(MpiCountSeik(1))
                MpiCountSeik=SeikDim
                MpiCountSeik(NotWorkingMember)=0
                
                allocate(MpiDisplacementSeik(0:SeikDim))
                MpiDisplacementSeik = huge(MpiDisplacementSeik(1))
                MpiDisplacementSeik(0)=0
                do i=1, SeikDim
                    MpiDisplacementSeik(i)=MpiDisplacementSeik(i-1)+MpiCountSeik(i-1)
                end do
            
            end if
            
            
            if(EnsembleRank==NotWorkingMember) then
                
                allocate(AllWeights(0:SeikDim))                    
                AllWeights = huge(AllWeights(0))
                AllWeights = 1.0d0/(SeikDim+1)
                
                allocate(TTTSeik(SeikDim,SeikDim))
                TTTSeik = huge(TTTSeik(1,1))
                call TTTSeik_builder()
                
                
                if(MyRank==0) then
            
                    allocate(CovSeik1(SeikDim,SeikDim))
                    CovSeik1 = huge(CovSeik1(1,1))
                    
                    allocate(CovSmoother1part(SeikDim,SeikDim))
                    CovSmoother1part = huge(CovSmoother1part(1,1))
                
                end if
                
            end if
            
      end subroutine 
      
      subroutine TTTSeik_builder()
            implicit none
            double precision, dimension(0:SeikDim,SeikDim) :: matrixT
            double precision, dimension(SeikDim,0:SeikDim) :: matrixT2
            integer :: indexi
            
            do indexi=1, SeikDim
                matrixT(:,indexi)=-AllWeights
                matrixT(indexi-1,indexi)=matrixT(indexi-1,indexi)+1
            end do
            do indexi=1, SeikDim
                TTTSeik(indexi,:)=matmul(AllWeights*matrixT(:,indexi),matrixT)
            end do
            
      end subroutine

      subroutine clean_seik()
            deallocate(trnEnsemble)
            deallocate(trnEnsembleWeighted)
            deallocate(BaseMember)
            deallocate(ModelErrorDiag)
            deallocate(LSeik)
            deallocate(MpiCount)
            deallocate(MpiDisplacement)
            deallocate(TempVecSeik)
            deallocate(TempSliceSeik)
            deallocate(TempMatrixSeik)
            
            if (MyRank==0) then
                deallocate(TempSliceSeik2)
                deallocate(MpiCountSeik)
                deallocate(MpiDisplacementSeik)
            end if
            
            if (EnsembleRank==0) then
                deallocate(AllWeights)
                deallocate(TTTSeik)
                
                if (MyRank==0) then
                    deallocate(CovSeik1)
                    deallocate(CovSmoother1part)
                end if
                
            end if
            
            
       
      end subroutine

      END MODULE memory_seik

