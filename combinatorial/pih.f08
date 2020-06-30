module util
    implicit none

    contains

    function arange(n) result(a)
    implicit none
    integer, intent(in) :: n
    integer :: i, a(n)
    do i = 1, n
        a(i) = i
    end do
    end function

    subroutine swap_int(a, b)
    implicit none
    integer, intent(inout) :: a, b
    integer :: c
    c = a
    a = b
    b = c
    end subroutine

    subroutine write_matrix(unit, a)
    implicit none
    integer, intent(in) :: unit, a(:, :)
    integer :: i
    character(len=7) :: fmt
    write(unit, '("shape",1x,2i6)') shape(a)
    write(fmt, '("(",i3,"i4)")') size(a, 2)
    do i = 1, size(a, 1)
        write(unit, fmt) a(i, :)
    end do
    end subroutine

end module 


module combinations
    implicit none

    contains

    logical &
    function lex_gt(a, b)
    implicit none
    integer, intent(in) :: a(:), b(:)
    integer :: i
    lex_gt = .false.
    do i = 1, size(a)
        if (a(i) < b(i)) then
            lex_gt = .true.
            return
        else if (a(i) > b(i)) then
            return
        end if
    end do
    end function

    subroutine symmetries(a, c)
    implicit none
    integer, intent(in) :: a(:, :)
    integer, intent(out) :: c
    integer :: s(size(a, 1)), o(size(a, 1)) 
    integer :: m, n, k, sk, unit

    m = size(a, 1)
    n = size(a, 1)
    o = 0
    s = 0
    k = 1
    c = 0
    open(newunit=unit, file='symmetries.txt')

100 sk = s(k) + 1
    if (sk > n) goto 102
    s(k) = sk
    if (o(sk) > 0) goto 100
    if (any(a(1:k-1, k) /= a(s(1:k-1), s(k)))) goto 100
    o(sk) = 1
    if (k < m) then 
        k = k + 1
        goto 100
    else
        c = c + 1
        write(unit, '(i12,6x,100i4)') c, s
    end if

101 o(s(k)) = 0
    goto 100

102 s(k) = 0
    k = k - 1
    if (k == 0) goto 103
    goto 101

103 close(unit)
    write(*, *) 'symmetries:', c
    end

    subroutine embeddings(t, a, isym)
    implicit none
    integer, intent(in) :: t(:), a(:, :)
    logical, intent(in), optional :: isym
    logical :: use_sym = .true.
    integer, allocatable :: sym(:, :)
    integer :: s(size(t)), o(size(a, 1)) 
    integer :: m, n, j, k, sk, c, nsym, unit

    m = size(t)
    n = size(a, 1)

    if (present(isym)) use_sym = isym

    if (use_sym) then
        call symmetries(a, nsym)
        allocate(sym(nsym, n))
        open(newunit=unit, file='symmetries.txt')
        do k = 1, nsym
            read(unit, *) j, sym(k, :)
        end do
        close(unit)
    end if

    o = 0
    s = 0
    k = 1
    c = 0
    open(newunit=unit, file='embeddings.txt')

100 sk = s(k) + 1
    if (sk > n) goto 102
    s(k) = sk
    if (o(sk) > 0) goto 100
    do j = 1, k - 1
       if (t(j) == t(k) .and. sk < s(j)) goto 100
    end do
    if (use_sym) then
        do j = 1, nsym
            if (lex_gt(sym(j, s(1:k)), s(1:k))) goto 100
        end do
    end if
    o(sk) = 1
    if (k < m) then 
        k = k + 1
        goto 100
    else
        c = c + 1
        write(unit, '(i12,6x,100i4)') c, s
    end if

101 o(s(k)) = 0
    goto 100

102 s(k) = 0
    k = k - 1
    if (k == 0) goto 103
    goto 101

103 close(unit)
    write(*, *) 'embeddings:', c
    end

end module


module combinations_tests
    use util
    use combinations
    implicit none

    contains

    subroutine chain(n, a)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: a(n, n)
    integer :: i
    a = 0
    do i = 1, n-1
        a(i, i+1) = 1
        a(i+1, i) = 1
    end do
    end

    subroutine ring(n, a)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: a(n, n)
    call chain(n, a)
    a(1, n) = 1
    a(n, 1) = 1
    end

    subroutine test_symmetries
    use util, only: write_matrix
    use combinations, only: embeddings
    implicit none
    integer, parameter :: n=5
    integer :: a(n, n)
    integer :: nsym
    call chain(n, a)
    call write_matrix(0, a)
    call symmetries(a, nsym)
    end 


    subroutine test_embeddings
    use combinations, only: embeddings
    implicit none
    integer, parameter :: n=4, m=2
    integer :: t(m), a(n, n), i
    t = arange(m)
    !t = 1
    call ring(n, a)
    call embeddings(t, a)
    end 

    subroutine test_all
    call test_symmetries
    call test_embeddings
    end subroutine

end module


    program main
    use combinations
    use util
    implicit none
    integer :: n, i, j, k, l, unit, a, t, io
    allocatable :: a(:, :), t(:), io(:)
    open(newunit=unit, file='graph.txt')
    read(unit, *) n
    allocate(a(n, n), t(n), io(n))
    a = 0
    do i = 1, n
        read(unit, *) j, t(i), k, (io(l), l=1, k)
        if (i /= j) then
            print*, 'input error!'
            stop
        end if
        do l = 1, k
            a(i, io(l)) = 1
            a(io(l), i) = 1
        end do
    end do
    close(unit)
    !call write_matrix(0, a)
    !print*, t
    call embeddings(t, a)
    end program
