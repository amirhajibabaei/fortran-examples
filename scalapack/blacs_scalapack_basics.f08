module dummy
    use, intrinsic:: iso_fortran_env, only: output_unit, input_unit, error_unit
    implicit none

    type process
        integer             :: rank
        integer             :: nprocs
        integer             :: row, col
        integer             :: nrows, ncols
        integer             :: ctxt
        contains
            procedure       :: print
    end type process

    contains


        subroutine print(p)
        implicit none
        class(process), intent(in) :: p
        write( output_unit, '("process",i3,2x,"out of",i3,2x)', advance='no' )  p%rank, p%nprocs
        write( output_unit, '("indexed",2i3,2x,"in process grid of the size",2i3)' )  p%row, p%col, p%nrows, p%ncols
        end subroutine

end module dummy



program test
use dummy
implicit none
type(process) :: proc

call blacs_pinfo( proc%rank, proc%nprocs )
call sl_init( proc%ctxt, 2, 2 )
call blacs_gridinfo( proc%ctxt, proc%nrows, proc%ncols, proc%row, proc%col )

call proc%print()

call blacs_gridexit( proc%ctxt )
call blacs_exit( 0 )
end program

