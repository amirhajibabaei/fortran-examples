program test
use easy_scalapack
implicit none
integer, parameter :: p = 2, q = 2, m = 9, n = 9, mb = 2, nb = 2
type(process_grid) :: pc
type(dense_matrix) :: dm, b
type(subarray)     :: sa, sb
integer            :: i, j, info
real               :: rnd

pc = process_grid(p,q)
dm = dense_matrix( m, n, mb, nb, pc )
b = dense_matrix( m, 1, mb, nb, pc )

!call dm%print_mapping()

do i = 1, dm%mg
    do j = 1, dm%ng
        call random_number(rnd)
        call dm%set( i, j, rnd )
    end do
    call b%set( i, 1, 1.0 )
end do

!call dm%print_matrix(6)

sa = subarray( dm )

! LU factorization
call psgetrf( sa%m, sa%n, sa%la, sa%ia, sa%ja, sa%desc, sa%ipiv, info )
if( info /= 0 ) call pc%stop_all("LU factorization unsuccessful")

sb = subarray(b)
call psgetrs( 'N', sa%m, 1, sa%la, sa%ia, sa%ja, sa%desc, sa%ipiv, sb%la, sb%ia, sb%ja, sb%desc, info )

!call psgetri( sa%n, sa%la, sa%ia, sa%ja, sa%desc, sa%ipiv, sa%work, sa%lwork, sa%iwork, sa%liwork, info ) still problematic!

call b%print_matrix(6)

call pc%destroy()
call pc%stop_all()

end program

