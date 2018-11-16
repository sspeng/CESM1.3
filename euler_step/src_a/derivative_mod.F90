module derivative_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np, nc, npdg, nep
implicit none

  type, public :: derivative_t
     real (kind=real_kind) :: Dvv(np,np)
     real (kind=real_kind) :: Dvv_diag(np,np)
     real (kind=real_kind) :: Dvv_twt(np,np)
     real (kind=real_kind) :: Mvv_twt(np,np)  ! diagonal matrix of GLL weights
     real (kind=real_kind) :: Mfvm(np,nc+1)
     real (kind=real_kind) :: Cfvm(np,nc)
     real (kind=real_kind) :: Sfvm(np,nep)
     real (kind=real_kind) :: legdg(np,np)
  end type derivative_t

end module derivative_mod
