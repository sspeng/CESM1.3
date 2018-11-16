module control_mod
  use kinds, only : real_kind

  real (kind=real_kind), public :: hypervis_power=0     ! if not 0, use variable hyperviscosity based on element area
  real (kind=real_kind), public :: hypervis_scaling=0      ! use tensor hyperviscosity



  integer, public, parameter :: west  = 1
  integer, public, parameter :: east  = 2
  integer, public, parameter :: south = 3
  integer, public, parameter :: north = 4

  integer, public, parameter :: swest = 5
  integer, public, parameter :: seast = 6
  integer, public, parameter :: nwest = 7
  integer, public, parameter :: neast = 8





end module control_mod
