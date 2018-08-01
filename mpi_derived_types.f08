!
! author: Amir H. T. 
!
module glob
      integer, parameter :: m=3, n=2, num_items=3

      type dummy_type
              integer    :: a(m)
              real       :: b(n)
              integer    :: c 
      end type dummy_type
end module glob


program test
      use mpi_f08
      use glob
      implicit none
      integer                   :: i, rank, ierr, tag, lengths(num_items)
      type(mpi_datatype)        :: my_mpi_type, datatypes(num_items)
      type(dummy_type)          :: foo
      integer(mpi_address_kind) :: locs(num_items)

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,rank,ierr)

      ! initiate dummy_type
      if( rank==0 ) then
              foo%a = [(i,i=1,m)]
              foo%b = [(real(i),i=1,n)]
              foo%c = 7
      elseif( rank==1 ) then
              foo%a = [(10+i,i=1,m)]
              foo%b = [(real(10+i),i=1,n)]
              foo%c = 17
      end if
      print*, "before:", rank, foo

      ! create an mpi_type for the dummy_type
      call mpi_get_address(foo%a,locs(1),ierr); lengths(1) = m; datatypes(1) =  mpi_integer
      call mpi_get_address(foo%b,locs(2),ierr); lengths(2) = n; datatypes(2) =  mpi_real
      call mpi_get_address(foo%c,locs(3),ierr); lengths(3) = 1; datatypes(3) =  mpi_integer
      locs = locs - locs(1)
      call mpi_type_create_struct(num_items,lengths,locs,datatypes,my_mpi_type,ierr)
      call mpi_type_commit(my_mpi_type,ierr)

      ! send-recv dummy_type
      tag = 12345
      if( rank==0 ) then
              call mpi_send(foo%a,1,my_mpi_type,1,tag,mpi_comm_world,ierr)
      elseif( rank==1 ) then
              call mpi_recv(foo%a,1,my_mpi_type,0,tag,mpi_comm_world,mpi_status_ignore,ierr)
      end if
      print*, " after:", rank, foo
  
      call mpi_finalize(ierr)
end program

