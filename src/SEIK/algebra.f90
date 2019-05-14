subroutine InvMatMul(A,B,nside,ierr)
    ! The upper triangular part of A will be substituted with the chiolesky decomposition of A, the result X of AX=B will be written over B
    implicit none
    integer, intent(in) :: nside
    integer, intent(out) :: ierr
    double precision, dimension (nside,nside), intent(inout) :: A, B 
    !double precision, dimension (nside,nside), intent(out) :: X
    
    !call dpotrf( 'U', nside, A, nside, ierr )
    !if (ierr.ne.0) error stop 'Matrix is numerically singular!'
    !call dpotrs( 'U', nside, nside, A, nside, B, nside, ierr)
    !if (ierr.ne.0) error stop 'Matrix inversion failed!'
    call dposv( 'U', nside, nside, A, nside, B, nside, ierr)
    if (ierr.ne.0) error stop 'Matrix inversion failed!'
end subroutine InvMatMul

subroutine OrtMatrix(CurrentMatrix, ColumnSize, RowSize, nRows)
    implicit none
    integer, intent(in) :: ColumnSize, RowSize, nRows
    double precision, dimension (:, :), intent(inout) :: CurrentMatrix
    integer :: indexi, indexj
    
    do indexi=1, nRows
        call NormalNumber(CurrentMatrix(:,RowSize+indexi),ColumnSize)
        CurrentMatrix(:,RowSize+indexi)=CurrentMatrix(:,RowSize+indexi)/norm2(CurrentMatrix(:,RowSize+indexi))
        do indexj=1, RowSize+indexi-1
            CurrentMatrix(:,RowSize+indexi)=CurrentMatrix(:,RowSize+indexi) & 
                -dot_product(CurrentMatrix(:,RowSize+indexi),CurrentMatrix(:,indexj))*CurrentMatrix(:,indexj)
            CurrentMatrix(:,RowSize+indexi)=CurrentMatrix(:,RowSize+indexi)/norm2(CurrentMatrix(:,RowSize+indexi))
        end do
    end do
end subroutine

subroutine NormalNumber(InputVector, nsize)
    implicit none
    integer, intent(in) :: nsize
    double precision, dimension(nsize), intent(out) :: InputVector
    double precision, dimension(2) :: workvector
    integer :: counter
    double precision :: circlerange
    
    counter=0
    do while (counter<nsize)
        call random_number(workvector)
        workvector=workvector*2.0d0
        workvector=workvector-1.0d0
        circlerange=workvector(1)*workvector(1)+workvector(2)*workvector(2)
        if ((circlerange >= 1.0d0).or.(circlerange==0.0d0)) cycle
        circlerange=sqrt(-2.0d0*log(circlerange)/circlerange)
        workvector=workvector*circlerange
        counter=counter+1
        InputVector(counter)=workvector(1)
        if (counter<nsize) then
            counter=counter+1
            InputVector(counter)=workvector(2)
        end if
    end do
end subroutine
