program test
use iso_fortran_env, only: real64, real32
!use f95_precision
use blas95
use lapack95
implicit none
integer, parameter :: n = 1000, pr = real64
real(pr)           :: a(n,n), b(n,n), c(n,n)
real(pr)           :: x(n), sigma 
integer            :: i, j, info, ipiv(n)

! construct a symmetric positive definite matrix
call random_number(x)
x = n*x
sigma = 1.0
do i = 1, n
   a(i,i) = 1.0
end do
do i = 1, n-1
   do j = i+1, n
      a(i,j) = exp( -(x(i)-x(j))**2/(2*sigma**2) )
      a(j,i) = a(i,j) 
   end do
end do

! b is going to be inverse of a
b = a

! use lapack for inversion 
call getrf(b,ipiv,info)
if( info/=0 ) call sysexit('LU decompositions not successful')
call getri(b,ipiv,info)
if( info/=0 ) call sysexit('Inversion not successful')

! c = I
!c = matmul( a, b )
call gemm( a, b, c )

! check if c = I
do i = 1, n
   c(i,i) = c(i,i) - 1.0
end do
write(*,*) 'residual:', sum( c**2 )/n**2

contains

        subroutine show_matrix(a)
        implicit none
        real a(:,:)
        integer i
        do i = 1, size(a,1)
           write(*,*) ( a(i,j), j=1,size(a,2) )
        end do
        end subroutine

        subroutine sysexit(input)
        class(*), intent(in) :: input
        select type(input)
        type is (character(len=*))
                print*, input
        end select
        stop
        end subroutine

end program
