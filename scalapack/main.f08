program test
use easy_scalapack
implicit none
integer, parameter :: p = 2, q = 2, m = 9, n = 9, mb = 2, nb = 2
type(process_grid) :: pc
type(dense_matrix) :: dm
integer            :: i, j, info
real               :: rnd

pc = process_grid(p,q)
dm = dense_matrix( m, n, mb, nb, pc )
call dm%print_mapping()
do i = 1, dm%mg
    do j = 1, dm%ng
        call random_number(rnd)
        call dm%set( i, j, rnd )
    end do
end do

!call dm%print_matrix(6)

! LU factorization
call psgetrf( dm%mg, dm%ng, dm%la, 1, 1, dm%desc, dm%ipiv, info )
if( info /= 0 ) call pc%stop_all("LU factorization unsuccessful")

! Inversion
! not completed yet! call psgetri( dm%mg, dm%la, 1, 1, dm%desc, work, lwork, dm%ipiv, info )


call pc%destroy()
call pc%stop_all()

end program

