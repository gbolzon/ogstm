module StringSEIK
implicit none
contains
function str2int(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: str2int, ierr

    read(str,*,iostat=ierr) str2int
    if (ierr/=0) error stop "not a valid string in str2int"
end function

function int2str(number, ndigits)
    implicit none

    integer, intent(in) :: number, ndigits
    character(len=ndigits) :: int2str
    integer :: ierr
    character(len=10) :: str

    write(str,"(I0)") ndigits
    write(int2str,"(I0."//trim(str)//")",iostat=ierr) number
    if (ierr/=0) error stop "not a valid number in int2str"
end function
end module
