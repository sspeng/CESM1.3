module hybvcoord_mod
use kinds, only: r8 => real_kind
use dimensions_mod, only: plev => nlev, plevp => nlevp
implicit none
type, public :: hvcoord_t
  real(r8) ps0          ! base state surface-pressure for level definitions
  real(r8) hyai(plevp)  ! ps0 component of hybrid coordinate - interfaces
  real(r8) hyam(plev)   ! ps0 component of hybrid coordinate - midpoints
  real(r8) hybi(plevp)  ! ps  component of hybrid coordinate - interfaces
  real(r8) hybm(plev)   ! ps  component of hybrid coordinate - midpoints
  real(r8) hybd(plev)   ! difference in b (hybi) across layers
  real(r8) prsfac       ! log pressure extrapolation factor (time, space independent)
  real(r8) etam(plev)   ! eta-levels at midpoints
  real(r8) etai(plevp)  ! eta-levels at interfaces
  integer  nprlev       ! number of pure pressure levels at top
  integer  pad
end type
end module hybvcoord_mod
