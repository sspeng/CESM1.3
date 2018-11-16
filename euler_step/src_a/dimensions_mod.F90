module dimensions_mod
  implicit none
  integer, parameter, public :: np = 4
  integer, parameter, public :: nlev = 30
  integer, parameter, public :: npdg = 0
  integer, parameter, public :: nc = 4
  integer, parameter, public :: nep = 9
  integer, parameter, public :: nelemd = 43
  integer, parameter, public :: qsize = 25
  integer, parameter, public :: nlevp = 31
  integer, parameter, public :: qsize_d = 25
  integer, public  :: max_elements_attached_to_node = 4
  integer, public  :: s_nv = 6
  integer, public  :: max_corner_elem               = 1 !max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8 !4 + 4*max_corner_elem
end module dimensions_mod
