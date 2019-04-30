      MODULE memory_seik
       
      USE modul_param

      IMPLICIT NONE

      public

      INTEGER :: SeikDim ! Reduced Dimension of the error subspace.
      INTEGER :: EnsembleComm, EnsembleRank, EnsembleSize
      double precision, allocatable, dimension (:,:,:,:) :: WeightSeik, trnEnsemble, trnEnsembleWeighted
       
      CONTAINS
       
      subroutine myalloc_seik()
            allocate(trnEnsemble(jpk,jpj,jpi,jptra))                    
            trnEnsemble = huge(trnEnsemble(1,1,1,1))
            
            allocate(WeightSeik(jpk,jpj,jpi,jptra))                    
            WeightSeik = huge(WeightSeik(1,1,1,1))
            WeightSeik = 1.0d0/(SeikDim+1)

            allocate(trnEnsembleWeighted(jpk,jpj,jpi,jptra))
            trnEnsembleWeighted = huge(trnEnsembleWeighted(1,1,1,1))
            
      end subroutine

      subroutine clean_seik()
            deallocate(trnEnsemble)
            deallocate(WeightSeik)
            deallocate(trnEnsembleWeighted)
       
      end subroutine

      END MODULE memory_seik

