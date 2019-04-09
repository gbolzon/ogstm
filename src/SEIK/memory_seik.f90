       MODULE memory_seik
       
       USE modul_param

       IMPLICIT NONE

       public

       INTEGER :: SeikDim ! Reduced Dimension of the error subspace.
       INTEGER :: EnsembleComm, EnsembleRank, EnsembleSize
       double precision, allocatable, dimension (:,:,:,:) :: WeightSeik, trnEnsemble
       
       CONTAINS
       
       subroutine myalloc_seik()
            allocate(trnEnsemble(jpk,jpj,jpi,jptra))                    
            trnEnsemble = huge(trnEnsemble(1,1,1,1))
            
            allocate(WeightSeik(jpk,jpj,jpi,jptra))                    
            WeightSeik = huge(WeightSeik(1,1,1,1))
            
       end subroutine
       

       END MODULE memory_seik

