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
            procedure       :: print, destroy, stop_all
    end type 

    interface process_grid
        procedure           :: init_processes
    end interface 


    type dense_matrix
        type(process_grid)  :: proc
        integer             :: descriptor(9)
        integer             :: rows_glob, cols_glob
        integer             :: rows_block, cols_block
        integer             :: rows_local, cols_local
        real, allocatable   :: localarray(:,:)
        contains
            procedure       :: bc_coords => block_cyclic_coords
            procedure       :: set       => set_dense_matrix
    end type 

    interface dense_matrix
        procedure           :: init_dense_matrix
    end interface


    contains

        ! ###########################################################
        ! ########################################################### process grid
        ! ###########################################################

        function init_processes( p, q ) result(proc)
        implicit none
        integer, intent(in) :: p, q
        type(process_grid)  :: proc
        external            :: blacs_pinfo, sl_init, blacs_gridinfo
        call blacs_pinfo( proc%rank, proc%nprocs )
        if( p*q/=proc%nprocs ) call proc%stop_all( "p*q /= nprocs!" )
        call sl_init( proc%ctxt, p, q )
        call blacs_gridinfo( proc%ctxt, proc%nrows, proc%ncols, proc%row, proc%col )
        proc%row_1st = 0
        proc%col_1st = 0
        end function

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine destroy( proc )
        implicit none
        class(process_grid), intent(in) :: proc
        external                        :: blacs_gridexit
        call blacs_gridexit( proc%ctxt )
        end subroutine

        subroutine stop_all( proc, message )
        implicit none
        class(process_grid), intent(in)        :: proc
        character(len=*), intent(in), optional :: message
        external                               :: blacs_exit
        if( present(message) .and. proc%rank==0 ) print *, message
        call blacs_exit(0)
        stop
        end subroutine

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine print(p)
        implicit none
        class(process_grid), intent(in) :: p
        write( output_unit, '("process",i3,2x,"out of",i3,2x)', advance='no' )  p%rank, p%nprocs
        write( output_unit, '("indexed",2i3,2x)', advance='no' )  p%row, p%col 
        write( output_unit, '("in process grid of the size",2i3)') p%nrows, p%ncols
        end subroutine

        ! ###########################################################
        ! ########################################################### dense matrix
        ! ###########################################################

        function init_dense_matrix( rows_glob, cols_glob, rows_block, cols_block, proc ) result(dm)
        implicit none
        integer, intent(in)             :: rows_glob, cols_glob
        integer, intent(in)             :: rows_block, cols_block
        class(process_grid), intent(in) :: proc
        type(dense_matrix)              :: dm
        integer                         :: info, numroc
        external                        :: numroc, descinit
        dm%rows_glob = rows_glob
        dm%cols_glob = cols_glob
        dm%rows_block = rows_block
        dm%cols_block = cols_block
        dm%rows_local = max( 1,  numroc( dm%rows_glob, dm%rows_block, proc%row, proc%row_1st, proc%nrows ) )
        dm%cols_local = max( 1,  numroc( dm%cols_glob, dm%cols_block, proc%col, proc%col_1st, proc%ncols ) )
        call descinit( dm%descriptor , dm%rows_glob, dm%cols_glob, dm%rows_block, dm%cols_block, proc%row_1st, proc%col_1st, proc%ctxt, dm%rows_local, info )
        if( info/=0 ) call proc%stop_all( "descinit info /=0" )
        dm%proc = proc
        allocate(dm%localarray(dm%rows_local,dm%cols_local))
        end function

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine block_cyclic_1d( iglob, nblock, ifirst, nproc, iproc, iloc )
        implicit none
        integer, intent(in)  :: iglob, nblock, ifirst, nproc
        integer, intent(out) :: iproc, iloc
        iproc = mod( ifirst + (iglob-1)/nblock , nproc )
        iloc = mod( iglob-1, nblock ) + 1
        end subroutine


        subroutine block_cyclic_coords( dm, iglob, jglob, iproc, jproc, iloc, jloc ) 
        implicit none
        class(dense_matrix), intent(in) :: dm
        integer,            intent(in)  :: iglob, jglob
        integer,            intent(out) :: iproc, jproc, iloc, jloc
        call block_cyclic_1d( iglob, dm%rows_block, dm%proc%row_1st, dm%proc%nrows, iproc, iloc)
        call block_cyclic_1d( jglob, dm%cols_block, dm%proc%col_1st, dm%proc%ncols, jproc, jloc)
        end subroutine

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        subroutine set_dense_matrix( dm, iglob, jglob, val )
        implicit none
        class(dense_matrix), intent(inout) :: dm
        integer,                intent(in) :: iglob, jglob
        real,                   intent(in) :: val
        integer                            :: iproc, jproc, iloc, jloc
        call dm%bc_coords( iglob, jglob, iproc, jproc, iloc, jloc )
        if( dm%proc%row==iproc .and. dm%proc%col==jproc ) then
            dm%localarray( iloc, jloc ) = val
        end if
        end subroutine

end module dummy



program test
use dummy
implicit none
integer, parameter :: p = 2, q = 2, m = 9, n = 9, mb = 2, nb = 2
type(process_grid) :: proc
type(dense_matrix) :: dm
integer            :: i, j, ii, jj, prow, pcol
integer, external           :: blacs_pnum

proc = process_grid(p,q)
dm = dense_matrix( m, n, mb, nb, proc )

do i = 1, m
   do j = 1, n
      call dm%bc_coords( i, j, prow, pcol, ii, jj )
      if(dm%proc%rank==0) write(*,*) i, j, blacs_pnum( dm%proc%ctxt, prow, pcol )
   end do
end do

call proc%destroy()
call proc%stop_all()

end program

