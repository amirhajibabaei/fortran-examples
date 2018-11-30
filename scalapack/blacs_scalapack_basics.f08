module dummy
    use, intrinsic:: iso_fortran_env, only: output_unit, input_unit, error_unit
    implicit none


    type process_grid
        integer             :: rank
        integer             :: nprocs
        integer             :: row, col
        integer             :: nrows, ncols
        integer             :: row_1st, col_1st
        integer             :: ctxt
        contains
            procedure       :: print, destroy
    end type 

    interface process_grid
        procedure           :: init_processes
    end interface 


    type dense_matrix
        integer             :: desca(9)
        integer             :: m, n, mb, nb
        integer             :: lld, lfd
    end type 

    interface dense_matrix
        procedure           :: init_dense_matrix
    end interface


    contains

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ process grid

        function init_processes( p, q ) result(proc)
        implicit none
        integer, intent(in) :: p, q
        type(process_grid)  :: proc
        external            :: blacs_pinfo, sl_init, blacs_gridinfo
        call blacs_pinfo( proc%rank, proc%nprocs )
        if( p*q/=proc%nprocs ) then
            print *, "p*q /= nprocs!"
            call blacs_exit(0)
            stop
        end if
        call sl_init( proc%ctxt, p, q )
        call blacs_gridinfo( proc%ctxt, proc%nrows, proc%ncols, proc%row, proc%col )
        proc%row_1st = 0
        proc%col_1st = 0
        end function

        subroutine destroy( proc )
        implicit none
        class(process_grid), intent(in) :: proc
        external                        :: blacs_gridexit
        call blacs_gridexit( proc%ctxt )
        end subroutine

        subroutine print(p)
        implicit none
        class(process_grid), intent(in) :: p
        write( output_unit, '("process",i3,2x,"out of",i3,2x)', advance='no' )  p%rank, p%nprocs
        write( output_unit, '("indexed",2i3,2x)', advance='no' )  p%row, p%col 
        write( output_unit, '("in process grid of the size",2i3)') p%nrows, p%ncols
        end subroutine

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ dense matrix

        function init_dense_matrix( m, n, proc ) result(dm)
        implicit none
        integer, intent(in)             :: m, n
        class(process_grid), intent(in) :: proc
        type(dense_matrix)              :: dm
        integer                         :: info, numroc
        external                        :: numroc, descinit
        dm%m = m
        dm%n = n
        dm%mb = 64
        dm%nb = 64
        dm%lld = max( 1,  numroc( dm%m, dm%mb, proc%row, proc%row_1st, proc%nrows ) )
        dm%lfd = max( 1,  numroc( dm%n, dm%nb, proc%col, proc%col_1st, proc%ncols ) )
        call descinit( dm%desca , dm%m, dm%n, dm%mb, dm%nb, proc%row_1st, proc%col_1st, proc%ctxt, dm%lld, info )
        if( info/=0 ) then 
            print*, "descinit info /=0"
        end if
        end function

end module dummy



program test
use dummy
implicit none
integer, parameter :: p = 2, q = 2, m = 1024, n = 1024
type(process_grid) :: proc
type(dense_matrix) :: dm
proc = process_grid(p,q)
call proc%print()
dm = dense_matrix( m, n, proc )
call proc%destroy()
call blacs_exit( 0 )
end program

