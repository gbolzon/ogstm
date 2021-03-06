       MODULE ZDF_mem 

       USE modul_param 
       USE myalloc

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public


      INTEGER :: dimen_jvzdf
      INTEGER, allocatable :: jarr_zdf(:,:),jarr_zdf_flx(:,:)
      double precision, allocatable :: zwd(:,:), zws(:,:), zwi(:,:)
      double precision, allocatable :: zwx(:,:), zwy(:,:), zwz(:,:), zwt(:,:)


!!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_ZDF()

#ifdef __OPENMP1
      INTEGER :: ntids, omp_get_max_threads
      EXTERNAL :: omp_get_max_threads
#else
      INTEGER :: ntids = 1
#endif
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
#endif
       dimen_jvzdf=0

       allocate(jarr_zdf(2, jpi*jpj))    
        jarr_zdf     = huge(jarr_zdf(1,1))
       allocate(jarr_zdf_flx(jpi*jpj,jpk)) 
        jarr_zdf_flx = huge(jarr_zdf_flx(1,1))
       allocate(zwd(jpk, ntids))           
        zwd          = huge(zwd(1,1))
       allocate(zws(jpk, ntids))           
        zws          = huge(zws(1,1))
       allocate(zwi(jpk, ntids))           
        zwi          = huge(zwi(1,1))
       allocate(zwx(jpk, ntids))           
        zwx          = huge(zwx(1,1))
       allocate(zwy(jpk, ntids))           
        zwy          = huge(zwy(1,1))
       allocate(zwz(jpk, ntids))           
        zwz          = huge(zwz(1,1))
       allocate(zwt(jpk, ntids))           
        zwt          = huge(zwt(1,1))

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif


      END subroutine myalloc_ZDF
      
      
      
      subroutine clean_memory_zdf()

          deallocate(jarr_zdf)
          deallocate(jarr_zdf_flx)
          deallocate(zwd)
          deallocate(zws)
          deallocate(zwi)
          deallocate(zwx)
          deallocate(zwy)
          deallocate(zwz)
          deallocate(zwt)

      end subroutine clean_memory_zdf



      END MODULE 
