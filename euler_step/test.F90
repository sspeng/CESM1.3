program comp
implicit none
  real (kind=8) :: a, a1, a2, a3, a4
  a = 0.5
  a1 = 0.1
  a2 = 2.0
  a3 = 8.3
  a4 = 1.9
  a = MAX(a, a1, a2, a3, a4)
  print *, a
end program comp

subroutine func(test)
  real :: test(16*2)
  integer :: i
  do i = 1, 16*2
    print *, test(i)
  enddo

end subroutine func
