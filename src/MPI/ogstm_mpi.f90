
MODULE ogstm_mpi_module

#ifdef ExecDA
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include "petsc/finclude/petscvec.h"
#else
#include <petsc/finclude/petscvecdef.h>
#endif
#endif

USE myalloc
USE mpi

#ifdef ExecDA
use DA_mem
use mpi_str, only: Var3DCommunicator
use petscvec, only: PETSC_COMM_WORLD, PETSC_NULL_CHARACTER
#endif

implicit NONE



contains
!! mpp routines
!!
!! mynode
!! mpplnk_my
!! mpprecv
!! mppsend
!! mppwait
!! mppsync
!! mppstop
!!
!!!---------------------------------------------------------------------
!!!
!!!                       routine mynode
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Massively parallel processors
!!!     Find processor unit
!!!
!!   Input :
!!   -----
!!      argument                :
!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------

SUBROUTINE mynode

      INTEGER :: ierr
#ifdef ExecDA
      PetscErrorCode :: stat
#endif

#ifdef key_mpp_mpi
      CALL mpi_comm_rank(MPI_COMM_WORLD,GlobalRank,ierr)
      CALL mpi_comm_size(MPI_COMM_WORLD,GlobalSize,ierr)
      
      
#ifdef ExecDA
      call parlec ! in order to read DA_Nprocs and SeikDim

! 3DVar Code

      if(GlobalRank .lt. DA_Nprocs) then
        call MPI_Comm_split(MPI_COMM_WORLD, DA_Nprocs, GlobalRank, Var3DCommunicator, ierr)
        
        PETSC_COMM_WORLD = Var3DCommunicator
        call PetscInitialize(PETSC_NULL_CHARACTER,stat)
        CHKERRQ(stat)
      else
        call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, GlobalRank, Var3DCommunicator, ierr)
      endif
      
!SEIK Code

      if(GlobalRank .lt. GlobalSize-mod(GlobalSize,SeikDim+1)) then

        !call MPI_Comm_split(MPI_COMM_WORLD, mod(GlobalRank,SeikDim+1), GlobalRank, LocalComm, ierr)      
        !changed order in the communcator for better parallelization in I/O
        call MPI_Comm_split(MPI_COMM_WORLD, GlobalRank/(GlobalSize/(SeikDim+1)), GlobalRank, LocalComm, ierr)
        CALL mpi_comm_rank(LocalComm,myrank,ierr)
        CALL mpi_comm_size(LocalComm,CommSize,ierr)
        
        call MPI_Comm_split(MPI_COMM_WORLD, myrank, GlobalRank, EnsembleComm, ierr)
        CALL mpi_comm_rank(EnsembleComm,EnsembleRank,ierr)
        CALL mpi_comm_size(EnsembleComm,EnsembleSize,ierr)
        
        !if (EnsembleRank==NotWorkingMember) then
        !    call MPI_Comm_split(EnsembleComm, MPI_UNDEFINED, EnsembleRank, BaseComm, ierr)
        !else
        !    call MPI_Comm_split(EnsembleComm, 0, EnsembleRank, BaseComm, ierr)
        !end if
        
        ! error check
        if (.false.) then !this error check was prepared before to change the construction of LocalComm and EnsembleComm
        if((myrank .ne. GlobalRank/(SeikDim+1)) .or. (CommSize .ne. GlobalSize/(SeikDim+1)) .or. (EnsembleRank .ne. mod(GlobalRank,SeikDim+1)) .or. (EnsembleSize .ne. SeikDim+1)) then
          write(*,*)'Unexpected value! myrank = ',myrank,', expected = ',GlobalRank/(SeikDim+1),'.'
          write(*,*)'CommSize = ',CommSize,', expected = ',GlobalSize/(SeikDim+1),'.'
          write(*,*)'EnsembleRank = ',EnsembleRank,', expected = ',mod(GlobalRank,SeikDim+1),'.'
          write(*,*)'EnsembleSize = ',EnsembleSize,', expected = ',SeikDim+1,'.'
          write(*,*)'Something wrong in ogstm_mpi.f90, I will stop.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          error stop
        end if
        end if
      else
        if(GlobalRank == GlobalSize-mod(GlobalSize,SeikDim+1)) then
          write(*,*)'The total number of processes (Global size = ',GlobalSize,') is not a multiple of the number of ensemble members (SeikDim + 1 = ',SeikDim+1,').'
          write(*,*)'This code is under development, I am unable to manage this exception at this moment and I will stop.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          error stop
        endif
      endif
        

#else
      LocalComm = MPI_COMM_WORLD
      myrank = GlobalRank
      CommSize = GlobalSize
#endif

#else
      GlobalRank = 0
      GlobalSize = 1
      CommSize = 1
      myrank = 0
#endif

END SUBROUTINE

     
!!!---------------------------------------------------------------------
!!!
!!!                       routine mpplnk_my
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!      Message passing manadgement
!!!
!!   Method :
!!   -------
!!       Use mppsend and mpprecv function for passing mask between
!!       processors following neighboring subdomains.
!!
!!   Input :
!!   -----
!!      argument
!!              ptab            : variable array
!!
!!      common
!!            /COMDOM/ : domain parameters
!!                    nbondi : mark for "east-west local boundary"
!!                    nbondj : mark for "north-south local boundary"
!!                    noea   : number for local neighboring processors
!!                    nowe   : number for local neighboring processors
!!                    noso   : number for local neighboring processors
!!                    nono   : number for local neighboring processors
!!
!!----------------------------------------------------------------------
 SUBROUTINE mpplnk_my(ptab)
      implicit none
      double precision, dimension(jpk,jpj,jpi), intent(inout) :: ptab
      
      call mpplnkSeik(ptab)

#ifdef key_mpp_mpi_I_dont_want_this_part

      INTEGER jk,jj,ji
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER reqs3, reqs4, reqr3, reqr4
      INTEGER jw, packsize

!!     trcadvparttime = MPI_WTIME()

!!
!!2. East and west directions exchange
!!------------------------------------



!!
!!2.2 Migrations
!!
!!
!               3     4
!               |     ^
!               |     |
!               v     |
!           ________________
!          |                |
!     1<-- |                | 1 <--
!     2--> |                | 2 -->
!          |________________|
!               3     4
!               |     ^
!               |     |
!               v     |

    packsize=jpk*jpj

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain

          CALL mppsend(2,ptab(:,:,jpi-1),packsize,noea,0,reqs1)
          CALL mpprecv(1,ptab(:,:,  jpi),packsize,reqr1)

      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1, ptab(:,:    ,2),packsize,nowe,0,reqs1)
          CALL mppsend(2, ptab(:,:,jpi-1),packsize,noea,0,reqs2)

          CALL mpprecv(1,ptab(:,:,jpi),packsize,reqr1)
          CALL mpprecv(2,ptab(:,:,  1),packsize,reqr2)

      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain

          CALL mppsend(1,ptab(:,:,2), packsize, nowe,0, reqs1)
          CALL mpprecv(2,ptab(:,:,1), packsize, reqr1)


      ENDIF


!!
!!
!!3. North and south directions
!!-----------------------------
!!
!!3.1 Read Dirichlet lateral conditions
!!


      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
         DO jw=1,NORTH_count_send
              ji = NORTHpoints_send(1,jw)
              jk = NORTHpoints_send(2,jw)
              tn_send(jw) = ptab(jk,jpj-1,ji)
         ENDDO
     ENDIF
     IF(nbondj.eq.0.or.nbondj.eq.1) THEN
         DO jw=1,SOUTH_count_send
             ji = SOUTHpoints_send(1,jw)
             jk = SOUTHpoints_send(2,jw)
             ts_send(jw) = ptab(jk,2,ji)
         ENDDO


      ENDIF!    PACK_LOOP4


!!
!!2.2 Migrations
!!
!!

      IF(nbondj.eq.-1) THEN ! We are at the south side of the domain
          CALL mppsend(4,tn_send,NORTH_count_send,nono,0,reqs4)
          CALL mpprecv(3,tn_recv,NORTH_count_recv,reqr3)
          CALL mppwait(reqs4)
          CALL mppwait(reqr3)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(4, tn_send,NORTH_count_send,nono,0,reqs4)
          CALL mppsend(3, ts_send,SOUTH_count_send,noso,0,reqs3)
          CALL mpprecv(3,tn_recv,NORTH_count_recv,reqr3)
          CALL mpprecv(4,ts_recv,SOUTH_count_recv,reqr4)

          CALL mppwait(reqs4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr3)
          CALL mppwait(reqr4)
      ELSE IF(nbondj.eq.1) THEN ! We are at the north side of the domain
          CALL mppsend(3,ts_send, SOUTH_count_send, noso,0, reqs3)
          CALL mpprecv(4,ts_recv, SOUTH_count_recv, reqr4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr4)
      ENDIF



!!
!!2.3 Write Dirichlet lateral conditions
!!

      IF(nbondj.eq.0.or.nbondj.eq.1) THEN ! All but south boundary, we received from south

         DO jw=1,SOUTH_count_recv
              ji = SOUTHpoints_recv(1,jw)
              jk = SOUTHpoints_recv(2,jw)
              ptab(jk,1,ji)= ts_recv(jw)
         ENDDO

      ENDIF

      IF(nbondj.eq.-1.or.nbondj.eq.0) THEN ! All but north boundary, we received from north

        DO jw=1,NORTH_count_recv
              ji = NORTHpoints_recv(1,jw)
              jk = NORTHpoints_recv(2,jw)
             ptab(jk,jpj,ji)= tn_recv(jw)
        ENDDO

      ENDIF ! PACK_LOOP5


!!!  East - West waits

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF

#endif

END SUBROUTINE

SUBROUTINE mpplnkSeik(ptab) !this subroutine should be rewritten using black/white chess-like algorithm. whites send first, blacks recive first.
    !use ogstm_mpi_module
    !use MPI
    !use myalloc

      implicit none
      double precision, dimension(jpk,jpj,jpi), intent(inout) :: ptab


#ifdef key_mpp_mpi

      INTEGER jk,jj,ji
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER reqs3, reqs4, reqr3, reqr4
      INTEGER packsize

!!     trcadvparttime = MPI_WTIME()

!!
!!2. East and west directions exchange
!!------------------------------------



!!
!!2.2 Migrations
!!
!!
!               3     4
!               |     ^
!               |     |
!               v     |
!           ________________
!          |                |
!     1<-- |                | 1 <--
!     2--> |                | 2 -->
!          |________________|
!               3     4
!               |     ^
!               |     |
!               v     |

    packsize=jpk*jpj

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain

          CALL mppsend(2,ptab(:,:,jpi-1),packsize,noea,0,reqs1)
          CALL mpprecv(1,ptab(:,:,  jpi),packsize,reqr1)

      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1, ptab(:,:    ,2),packsize,nowe,0,reqs1)
          CALL mppsend(2, ptab(:,:,jpi-1),packsize,noea,0,reqs2)

          CALL mpprecv(1,ptab(:,:,jpi),packsize,reqr1)
          CALL mpprecv(2,ptab(:,:,  1),packsize,reqr2)

      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain

          CALL mppsend(1,ptab(:,:,2), packsize, nowe,0, reqs1)
          CALL mpprecv(2,ptab(:,:,1), packsize, reqr1)


      ENDIF

!!!  East - West waits

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF



!!
!!2.2 Migrations
!!
!!
      packsize=jpk*jpi
        BufferMPILinkSend3=ptab(:,jpj-1,:)        
        BufferMPILinkSend4=ptab(:,2,:)          
        ptab(:,1,:)=BufferMPILinkRecieve3                
        ptab(:,jpj,:)=BufferMPILinkRecieve4
        
      IF(nbondj.eq.-1) THEN ! We are at the south side of the domain
            BufferMPILinkSend4=ptab(:,2,:)    
          CALL mppsend(4,BufferMPILinkSend4,packsize,nono,0,reqs4)
          CALL mpprecv(3,BufferMPILinkRecieve3,packsize,reqr3)
          CALL mppwait(reqs4)
          CALL mppwait(reqr3)
            ptab(:,1,:)=BufferMPILinkRecieve3 
      ELSE IF(nbondj.eq.0) THEN
            BufferMPILinkSend3=ptab(:,jpj-1,:)        
            BufferMPILinkSend4=ptab(:,2,:)   
          CALL mppsend(4, BufferMPILinkSend4,packsize,nono,0,reqs4)
          CALL mppsend(3, BufferMPILinkSend3,packsize,noso,0,reqs3)
          CALL mpprecv(3,BufferMPILinkRecieve3,packsize,reqr3)
          CALL mpprecv(4,BufferMPILinkRecieve4,packsize,reqr4)

          CALL mppwait(reqs4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr3)
          CALL mppwait(reqr4)
            ptab(:,1,:)=BufferMPILinkRecieve3                
            ptab(:,jpj,:)=BufferMPILinkRecieve4
      ELSE IF(nbondj.eq.1) THEN ! We are at the north side of the domain
            BufferMPILinkSend3=ptab(:,jpj-1,:) 
          CALL mppsend(3,BufferMPILinkSend3,packsize, noso,0, reqs3)
          CALL mpprecv(4,BufferMPILinkRecieve4,packsize, reqr4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr4)
            ptab(:,jpj,:)=BufferMPILinkRecieve4
      ENDIF


#endif

END SUBROUTINE





!!!---------------------------------------------------------------------
!!!
!!!                       routine mppsend
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Send messag passing array
!!
!!   Input :
!!   -----
!!      argument                :
!!                   ktyp   -> Tag of the message
!!                   pmess  -> array of double precision to send
!!                   kbytes -> size of pmess in double precision
!!                   kdest  -> receive process number
!!                   kid    _> ? (note used)
!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------

SUBROUTINE mppsend(ktyp,pmess,kbytes,kdest,kid,ireqsend)


      double precision pmess(*)
      INTEGER kbytes,kdest,ktyp,kid, ireqsend

#ifdef key_mpp_mpi





      INTEGER iflag
      CALL mpi_isend(pmess,kbytes,mpi_real8,kdest,ktyp, &
     &    LocalComm,ireqsend,iflag)


#endif
      RETURN

END SUBROUTINE

!!!---------------------------------------------------------------------
!!!
!!!                       routine mpprecv
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Receive messag passing array
!!
!!   Input :
!!   -----
!!      argument
!!                   ktyp    -> Tag of the recevied message
!!                   pmess   -> array of double precision
!!                   kbytes  -> suze of the array pmess


!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------

 SUBROUTINE mpprecv(ktyp,pmess,kbytes,ireqrecv)

      double precision pmess(*)
      INTEGER   kbytes,ktyp, ireqrecv

#ifdef key_mpp_mpi

      INTEGER iflag

      CALL mpi_irecv(pmess,kbytes,mpi_real8,mpi_any_source,ktyp,LocalComm,ireqrecv,iflag)
#endif

      RETURN

END SUBROUTINE

!!!---------------------------------------------------------------------
!!!
!!!                       routine mppwait
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Wait message passing isend/irecv
!!
!!   Input :
!!   -----
!!      argument
!!----------------------------------------------------------------------
SUBROUTINE mppwait(req)
      INTEGER istatus(mpi_status_size), ierr
      integer req
      


#ifdef key_mpp_mpi
      call MPI_WAIT(req, istatus, ierr)
#endif
      RETURN
END SUBROUTINE

!!!---------------------------------------------------------------------
!!!
!!!                       routine mppsync
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Massively parallel processors, synchroneous
!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------
SUBROUTINE mppsync()

!!----------------------------------------------------------------------

#ifdef key_mpp_mpi

      INTEGER ierror

      CALL mpi_barrier(LocalComm,ierror) !sostituito. Il commento precedente era, per il momento lascio mpi_comm_world, da riverificare se si puo'/conviene sostituire con LocalComm

#endif
      RETURN
END SUBROUTINE

SUBROUTINE mppstop
!!!---------------------------------------------------------------------
!!!
!!!                       routine mppstop
!!!                     *******************
!!!
!!!  purpose :
!!!  --------
!!!     Stop massilively parallel processors method
!!


      INTEGER :: info
      INTEGER :: ierr

#ifdef key_mpp_mpi
      CALL mppsync
#ifdef ExecDA
! seik deallocation
      call mpi_comm_free(LocalComm, ierr)
      call mpi_comm_free(EnsembleComm, ierr)
#endif      
#endif

      RETURN
END SUBROUTINE

END MODULE ogstm_mpi_module
