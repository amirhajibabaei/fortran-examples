module dummy
    use, intrinsic:: iso_fortran_env, only: output_unit, input_unit, error_unit
    implicit none


    type process_grid
        integer              :: rank
        integer              :: size
        integer              :: myr, myc
        integer              :: npr, npc
        integer              :: r1st, c1st
        integer              :: ctxt
        contains
            procedure        :: print_process_grid, destroy, stop_all
    end type 

    interface process_grid
        procedure            :: init_processes
    end interface 


    type dense_matrix
        type(process_grid)   :: pc
        integer              :: desc(9)
        integer              :: mg, ng
        integer              :: mb, nb
        integer              :: ml, nl
        real, allocatable    :: la(:,:)
        integer, allocatable :: ipiv(:)
        contains
            procedure        :: locr, locc
            procedure        :: bc_coords => block_cyclic_coords
            procedure        :: set       => set_dense_matrix
            procedure        :: print_mapping, serial_matrix, print_matrix
    end type 

    interface dense_matrix
        procedure            :: init_dense_matrix
    end interface


    interface nums_lcm
        procedure least_common_multiple
    end interface

    interface nums_gcd
        procedure greatest_common_divisor
    end interface


    contains

        ! ###########################################################
        ! ########################################################### auxiliary 
        ! ###########################################################


        pure elemental &
        function least_common_multiple(x,y) result(c)
        implicit none
        integer, intent(in) :: x,y
        integer             :: c
        c = x*y / greatest_common_divisor(x,y)
        end function

        pure elemental &
        function greatest_common_divisor(x,y) result(c)
        implicit none
        integer, intent(in) :: x,y
        integer             :: a,b
        integer             :: t,c
        a = x
        b = y
        do while (b/=0)
            t = b
            b = mod(a,b)
            a = t
        end do
        c = abs(a)
        end function


        ! ###########################################################
        ! ########################################################### process grid
        ! ###########################################################

        function init_processes( p, q ) result(pc)
        implicit none
        integer, intent(in) :: p, q
        type(process_grid)  :: pc
        external            :: blacs_pinfo, sl_init, blacs_gridinfo
        call blacs_pinfo( pc%rank, pc%size )
        if( p*q/=pc%size ) call pc%stop_all( "p*q /= size!" )
        call sl_init( pc%ctxt, p, q )
        call blacs_gridinfo( pc%ctxt, pc%npr, pc%npc, pc%myr, pc%myc )
        pc%r1st = 0
        pc%c1st = 0
        end function

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine destroy( pc )
        implicit none
        class(process_grid), intent(in) :: pc
        external                        :: blacs_gridexit
        call blacs_gridexit( pc%ctxt )
        end subroutine

        subroutine stop_all( pc, message )
        implicit none
        class(process_grid), intent(in)        :: pc
        character(len=*), intent(in), optional :: message
        external                               :: blacs_exit
        if( present(message) .and. pc%rank==0 ) print *, message
        call blacs_exit(0)
        stop
        end subroutine

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine print_process_grid(p)
        implicit none
        class(process_grid), intent(in) :: p
        write( output_unit, '("process",i3,2x,"out of",i3,2x)', advance='no' )  p%rank, p%size
        write( output_unit, '("indexed",2i3,2x)', advance='no' )  p%myr, p%myc 
        write( output_unit, '("in process grid of the size",2i3)') p%npr, p%npc
        end subroutine

        ! ###########################################################
        ! ########################################################### dense matrix
        ! ###########################################################

        function init_dense_matrix( mg, ng, mb, nb, pc ) result(dm)
        implicit none
        integer, intent(in)             :: mg, ng
        integer, intent(in)             :: mb, nb
        class(process_grid), intent(in) :: pc
        type(dense_matrix)              :: dm
        integer                         :: info, numroc
        external                        :: numroc, descinit
        dm%mg = mg; dm%mb = mb
        dm%ng = ng; dm%nb = nb
        dm%ml = max( 1,  numroc( dm%mg, dm%mb, pc%myr, pc%r1st, pc%npr ) )
        dm%nl = max( 1,  numroc( dm%ng, dm%nb, pc%myc, pc%c1st, pc%npc ) )
        call descinit( dm%desc , dm%mg, dm%ng, dm%mb, dm%nb, pc%r1st, pc%c1st, pc%ctxt, dm%ml, info )
        dm%pc = pc
        if( info/=0 ) call pc%stop_all( "descinit info /=0" )
        allocate(dm%la(dm%ml,dm%nl))
        allocate(dm%ipiv(dm%ml+dm%mb))
        end function

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOCr(), LOCc()

        integer &
        function locr( dm, s ) 
        implicit none
        class(dense_matrix), intent(in) :: dm
        integer, intent(in), optional   :: s
        integer, external               :: numroc
        if( present(s) ) then
            locr = max( 1, numroc( s,     dm%mb, dm%pc%myr, dm%pc%r1st, dm%pc%npr ) )
        else
            locr = max( 1, numroc( dm%mg, dm%mb, dm%pc%myr, dm%pc%r1st, dm%pc%npr ) )
        end if
        end function

        integer &
        function locc( dm, s ) 
        implicit none
        class(dense_matrix), intent(in) :: dm
        integer, intent(in), optional   :: s
        integer, external               :: numroc
        if( present(s) ) then
            locc = max( 1,  numroc( s,     dm%nb, dm%pc%myc, dm%pc%c1st, dm%pc%npc ) )
        else
            locc = max( 1,  numroc( dm%ng, dm%nb, dm%pc%myc, dm%pc%c1st, dm%pc%npc ) )
        end if
        end function

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ i,j -> ip,jp and il,jl

        subroutine block_cyclic_1d( i, nb, p1st, np, ip, il )
        implicit none
        integer, intent(in)  :: i, nb, p1st, np
        integer, intent(out) :: ip, il
        ip = mod( p1st + (i-1)/nb , np )
        il = ((i-1)/(np*nb))*nb + mod( i-1, nb ) + 1
        end subroutine

        subroutine block_cyclic_coords( dm, i, j, ip, jp, il, jl ) 
        implicit none
        class(dense_matrix), intent(in) :: dm
        integer,            intent(in)  :: i, j
        integer,            intent(out) :: ip, jp, il, jl
        call block_cyclic_1d( i, dm%mb, dm%pc%r1st, dm%pc%npr, ip, il)
        call block_cyclic_1d( j, dm%nb, dm%pc%c1st, dm%pc%npc, jp, jl)
        end subroutine

        subroutine print_mapping(dm)
        implicit none
        class(dense_matrix), intent(in) :: dm
        integer                         :: i, j, ip, jp, il, jl, u
        if( dm%pc%rank==0 ) then
            open( newunit=u, file='dense_matrix_mapping.txt' )
            write(u,*) '#global', dm%mg, dm%ng
            write(u,*) '#block', dm%mb, dm%nb
            write(u,*) '#process', dm%pc%npr, dm%pc%npc
            do i = 1, dm%mg
                do j = 1, dm%ng
                    call dm%bc_coords(i,j,ip,jp,il,jl)
                    write(u,*) i,j,ip,jp,il,jl
                end do
            end do
            close(u)
        end if
        end subroutine

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine serial_matrix(dm)
        implicit none
        class(dense_matrix), intent(inout) :: dm
        integer                            :: i, j
        do i = 1, dm%mg
            do j = 1, dm%ng
                call dm%set(i,j, real(j-1)*dm%mg + real(i) )
            end do
        end do
        end subroutine


        subroutine print_matrix(dm,u)
        implicit none
        class(dense_matrix), intent(inout) :: dm
        integer,                intent(in) :: u
        real                               :: work( dm%mb )
        call pslaprnt( dm%mg, dm%ng, dm%la(1,1), 1, 1, dm%desc, 0, 0, "", u, work )
        ! note: pslaprnt for double, etc.
        end subroutine

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine set_dense_matrix( dm, i, j, val )
        implicit none
        class(dense_matrix), intent(inout) :: dm
        integer,                intent(in) :: i, j
        real,                   intent(in) :: val
        integer                            :: ip, jp, il, jl
        call dm%bc_coords( i, j, ip, jp, il, jl )
        if( dm%pc%myr==ip .and. dm%pc%myc==jp ) then
            dm%la( il, jl ) = val
        end if
        end subroutine

end module dummy



program test
use dummy
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

