program comp
implicit none
  !integer :: a = 5
  !integer :: b = 6
  !integer :: m = 7
  !integer :: c, d
  !c = min(a, b, m)
  !d = max(a, b, m)
  !print *, c, d
  !integer :: i, j, k
  !integer :: a(10, 10)
  !do j = 1, 10
  !  do i = 1, 10
  !    a(i, j) = j * 10 + i
  !  enddo
  !enddo
!
  !do j = 1, 10
  !  do i = 1, 10
  !    if ((i>3) .and. (i<5)) CYCLE
  !    print *, i, j, a(i, j)
  !    if (i>6) exit
  !  enddo
  !enddo
  !print *, 0D0
  !integer :: a(16), b(16), c(16)
  !integer :: i, j, k, sum_a = 0
  !do i = 1, 16
  !  a(i) = 3
  !  b(i) = 5
  !enddo
  !c = a*b
  !!do i = 1, 16
  !!  print *, c(i)
  !!enddo
  !sum_a = sum(a*b)
  !integer :: i, j, k
  !real :: mass
  !real :: cc(16), xx(16)
  !real :: test(4, 4, 2)
  !do i = 1, 16
  !  cc(i) = 5
  !  xx(i) = 5
  !enddo
  !mass = sum(cc*xx)
  !print *, mass
  !do k = 1, 2
  !  do j = 1, 4
  !    do i = 1, 4
  !      test(i, j, k) = 1
  !    enddo
  !  enddo
  !enddo
  !call func(test)
  !real :: a, b
  !a = -3.0d0
  !b = abs(a)
  !print *, a, b
  !real(kind=8) :: a, b
  !a = a + 1
  !print *, a
  integer :: i, j, k
  integer :: a(4, 4, 8), b(4, 4, 8), c(4,4)
  do k = 1, 8
    do j = 1, 4
      do i = 1, 4
        a(i,j,k) = 1
        b(i,j,k) = 1
      enddo
    enddo
  enddo
  do j = 1, 4
    do i = 1, 4
      c(i,j) = sum(a(i,j,:)*b(i,j,:))
    enddo
  enddo
  do j = 1, 4
    do i = 1, 4
      print *, c(i, j)
    enddo
  enddo
end program comp

subroutine func(test)
  real :: test(16*2)
  integer :: i
  do i = 1, 16*2
    print *, test(i)
  enddo

end subroutine func
