SUBROUTINE mpplnkSeik(ptab) !this subroutine should be rewritten using black/white chess-like algorithm. whites send first, blacks recive first.
    use ogstm_mpi_module
    use MPI
    use myalloc

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
      IF(nbondj.eq.-1) THEN ! We are at the south side of the domain
          CALL mppsend(4,ptab(:,2,:),packsize,nono,0,reqs4)
          CALL mpprecv(3,ptab(:,1,:),packsize,reqr3)
          CALL mppwait(reqs4)
          CALL mppwait(reqr3)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(4, ptab(:,2,:),packsize,nono,0,reqs4)
          CALL mppsend(3, ptab(:,jpj-1,:),packsize,noso,0,reqs3)
          CALL mpprecv(3,ptab(:,1,:),packsize,reqr3)
          CALL mpprecv(4,ptab(:,jpj,:),packsize,reqr4)

          CALL mppwait(reqs4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr3)
          CALL mppwait(reqr4)
      ELSE IF(nbondj.eq.1) THEN ! We are at the north side of the domain
          CALL mppsend(3,ptab(:,jpj-1,:),packsize, noso,0, reqs3)
          CALL mpprecv(4,ptab(:,jpj,:),packsize, reqr4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr4)
      ENDIF


#endif

END SUBROUTINE
