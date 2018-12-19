program main
use kinds , only: real_kind
use element_mod, only: element_t
use dimensions_mod, only: np, npdg, nlev, qsize
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use hybvcoord_mod, only:hvcoord_t
implicit none
  integer :: np1_qdp = 1
  integer :: n0_qdp = 2
  real(kind=real_kind) :: dt = 1
  type(element_t) :: elem(43)
  type(element_t) :: elem_test(43)
  type(hybrid_t) :: hybrid
  type(derivative_t) :: deriv
  type(hvcoord_t) :: hvcoord
  integer, parameter :: nets = 1
  integer, parameter :: nete = 43
  integer :: DSSopt = 3
  integer :: rhs_multiplier = 0

  real(kind=real_kind), dimension(np, np                        ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np, np                        ) :: div
  real(kind=real_kind), dimension(np, np, nlev                  ) :: dpdissk
  real(kind=real_kind), dimension(np, np, 2                     ) :: gradQ
  real(kind=real_kind), dimension(np, np, 2, nlev               ) :: Vstar
  real(kind=real_kind), dimension(np, np, nlev                  ) :: Qtens
  real(kind=real_kind), dimension(np, np, nlev                  ) :: dp, dp_star
  real(kind=real_kind), dimension(np, np, nlev, qsize, nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)                 :: DSSvar
  real(kind=real_kind) :: dp0(nlev), qim_val(nlev), qmax_val(nlev)
  integer :: ie, q, i, j, k, d, l
  integer :: rhs_viss = 0
  integer :: qbeg, qend, kbeg, kend
  integer :: kptr
  integer :: count
  external :: counter
  interface
    subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
      use kinds , only: real_kind
      use element_mod, only: element_t
      use hybvcoord_mod, only:hvcoord_t
      use hybrid_mod, only: hybrid_t
      use derivative_mod, only: derivative_t
      integer :: np1_qpd, n0_qdp, nets, nete, DSSopt, rhs_multiplier
      real(kind=real_kind) :: dt
      type(element_t) :: elem(:)
      type(hybrid_t) :: hybrid
      type(derivative_t) :: deriv
      type(hvcoord_t) :: hvcoord
    end subroutine euler_step
  end interface

  call athread_init()

  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,n0_qdp) = 1!ie * 1000000 + k * 100 + j * 10 + i + 0.2
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = 1!ie * 1000000 + k * 100 + j * 10 + i + 0.1
          enddo
        enddo
      enddo
    enddo
  enddo

  do ie = nets, nete
    do k = 1, nlev
      do j = 1, np
        do i = 1, np
          elem(ie)%derived%dp(i,j,k) = 1 !(ie * 1000000 + k * 100 + j * 10 + i)*2
          elem(ie)%derived%dpdiss_ave(i,j,k) = 1 !(ie * 1000000 + k * 100 + j * 10 + i)*2
          elem(ie)%derived%divdp_proj(i,j,k) = 1 ! ie * 1000000 + k * 100 + j * 10 + i
          elem(ie)%derived%divdp(i,j,k) = 1 !ie * 1000000 + k * 100 + j * 10 + i
          elem(ie)%derived%dpdiss_biharmonic(i,j,k) = 1! ie * 1000000 + k * 100 + j * 10 + i
        enddo
      enddo
    enddo
  enddo

  do ie = nets, nete
    do k = 1, nlev
      do j = 1, np
        do i = 1, np
          elem(ie)%derived%vn0(i,j,1,k) = 1!(ie * 1000000 + k * 100 + j * 10 + i)
          elem(ie)%derived%vn0(i,j,2,k) = 1!(ie * 1000000 + k * 100 + j * 10 + i)*2
        enddo
      enddo
    enddo
  enddo

  do l = 1, np
    do i = 1, np
      deriv%Dvv(i, l) = 10*l + i
    enddo
  enddo

  do ie = nets, nete
    do j = 1, np
      do i = 1, np
        elem(ie)%Dinv(i, j, 1, 1) = 1
        elem(ie)%Dinv(i, j, 2, 1) = 2
        elem(ie)%Dinv(i, j, 1, 2) = 3
        elem(ie)%Dinv(i, j, 2, 2) = 4
      enddo
    enddo
  enddo

  do ie = nets, nete
    do j = 1, np
      do i = 1, np
        elem(ie)%metdet(i, j) = 1
        elem(ie)%rmetdet(i, j) = 1
        elem(ie)%spheremp(i, j) = 1
      enddo
    enddo
  enddo

  do ie = nets, nete
    do j = 1, np
      do i = 1, np
        elem(ie)%variable_hyperviscosity(i,j) = 1
        elem(ie)%tensorVisc(i,j,1,1) = 1
        elem(ie)%tensorVisc(i,j,1,2) = 1
        elem(ie)%tensorVisc(i,j,1,2) = 1
        elem(ie)%tensorVisc(i,j,2,2) = 1
      enddo
    enddo
  enddo

  call euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )

  do i = 1, 5
    call counter(count)
  enddo
  !print *, count
end program main

subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
use kinds , only: real_kind
use physical_constants, only: rrearth
use element_mod, only: element_t
use dimensions_mod, only: np, npdg, nlev, qsize, qsize_d, max_corner_elem, max_neigh_edges, nelemd
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use edgetype_mod, only       : EdgeBuffer_t
implicit none
  integer, intent(in) :: np1_qdp, n0_qdp
  real(kind=real_kind) :: dt
  type(element_t), intent(inout) :: elem(:)
  type(hybrid_t), intent(in) :: hybrid
  type(derivative_t), intent(in) :: deriv
  type(hvcoord_t), intent(in) :: hvcoord
  integer, intent(in):: nets
  integer, intent(in) :: nete
  integer :: DSSopt
  integer :: rhs_multiplier
  real(kind=real_kind), dimension(np, np                           ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np, np                           ) :: div
  real(kind=real_kind), dimension(np, np,    nlev                  ) :: dpdissk
  real(kind=real_kind), dimension(np, np, 2                        ) :: gradQ
  real(kind=real_kind), dimension(np, np, 2, nlev                  ) :: Vstar
  real(kind=real_kind), dimension(np, np,    nlev                  ) :: Qtens
  real(kind=real_kind), dimension(np, np,    nlev                  ) :: dp, dp_star
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), dimension(np, np,    nlev,        nets:nete) :: dp_temp
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: dp_star_temp
  real(kind=real_kind), dimension(np, np, 2, nlev, qsize, nets:nete) :: gradQ_temp
  real(kind=real_kind), dimension(           nlev, qsize, nets:nete) :: qmax
  real(kind=real_kind), dimension(           nlev, qsize, nets:nete) :: qmin
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: Qtens_temp
  real(kind=real_kind), dimension(np, np,    nlev, qsize, nets:nete) :: Qtens_temp_test
  real(kind=real_kind), pointer, dimension(:,:,:) :: DSSvar
  real(kind=real_kind) :: dp0(nlev), qim_val(nlev), qmax_val(nlev)
  type (EdgeBuffer_t) :: edgeAdv, edgeAdvp1, edgeAdvQminmax, edgeAdv1,  edgeveloc
  integer :: ie, q, i, j, k, d
  integer :: limiter_option = 8
  integer :: rhs_viss = 1
  integer :: rkstage = 3
  integer :: qbeg, qend, kbeg, kend
  integer :: kptr
  real(kind=real_kind) :: time_s, time_e
  real(kind=real_kind) :: nu_p = 0.0D5
  real(kind=real_kind) :: nu_q = -1
  type(element_t) :: elem_test(43)
  integer :: count = 0
  real (kind=real_kind) :: tol_limiter = 5D-14

  !external :: slave_qdp_time_avg
  !type param_t
  !  integer*8 :: qdp
  !  integer :: rkstage, n0_qdp, np1_qdp, nets, nete, qsize \
  !     , qsize_d, step_elem
  !end type param_t
  !type(param_t) :: param_s
  !external :: slave_euler_step
  external :: slave_euler_step
  type param_t
    integer*8 :: qdp_s_ptr, qdp_leap_ptr,dp_s_ptr, divdp_proj_s_ptr   &
        , Qtens_biharmonic, qmax, qmin
    real(kind=real_kind) :: dt
    integer :: nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize, qsize_d
  end type param_t
  type(param_t) :: param_s

  external :: slave_euler_v
  type param_2d_t
    integer*8 :: qdp_s_ptr, qdp_leap_ptr, divdp_proj, dp, vn0, Dvv, Dinv         \
    , metdet, rmetdet, Qtens_biharmonic, divdp, dpdiss_biharmonic, spheremp      \
    , qmax, qmin
    real(kind=real_kind) :: dt, rrearth, nu_p, nu_q
    integer :: nets, nete, rhs_multiplier, qsize, qsize_d, n0_qdp, np1_qdp, limiter_option\
        , rhs_viss
  end type param_2d_t
  type(param_2d_t) :: param_2d_s

  external :: slave_compute_overlap
  type param_compute_t
    integer*8 :: Qtens_biharmonic, spheremp, dp0
    real(kind=real_kind) :: dt, nu_q
    integer :: nets, nete, qsize, step_elem, rhs_viss
  end type param_compute_t
  type(param_compute_t) :: param_compute_s

  external :: slave_rhs_multiplier
  type param_rsh_t
    integer*8 :: Qtens_biharmonic, dpdiss_ave, dp0
    integer :: nets, nete, qsize, step_elem
  end type param_rsh_t
  type(param_rsh_t) :: param_rhs_s

  interface
    subroutine biharmonic_wk_scalar(elem, Qtens_biharmonic, deriv, edgeAdv, hybrid, nets, nete)
      use kinds , only: real_kind
      use element_mod, only: element_t
      use hybrid_mod, only: hybrid_t
      use derivative_mod, only: derivative_t
      use edgetype_mod, only       : EdgeBuffer_t
      type(element_t) :: elem(:)
      real(kind=real_kind) :: Qtens_biharmonic(:,:,:,:,:)
      type(EdgeBuffer_t) :: edgeAdv
      type(derivative_t) :: deriv
      type(hybrid_t) :: hybrid
      integer :: nets, nete
    end subroutine biharmonic_wk_scalar
  end interface
interface
    subroutine neighbor_minmax(hybrid, edgeMinMax, nets, nete, min_neigh, max_neigh)
      use kinds , only: real_kind
      use hybrid_mod, only: hybrid_t
      use edgetype_mod, only       : EdgeBuffer_t
      real(kind=real_kind) :: min_neigh(:,:,:)
      real(kind=real_kind) :: max_neigh(:,:,:)
      type(EdgeBuffer_t) :: edgeMinMax
      type(hybrid_t) :: hybrid
      integer :: nets, nete
    end subroutine neighbor_minmax
  end interface

  interface
    subroutine edgeVpack(elem, edgeAdvp1, nets, nete, np1_qdp, DSSopt)
      use kinds , only: real_kind
      use element_mod, only: element_t
      use edgetype_mod, only: EdgeBuffer_t
      type(element_t) :: elem(:)
      type(EdgeBuffer_t) :: edgeAdvp1
      integer :: nets, nete, np1_qdp, DSSopt
    end subroutine edgeVpack
  end interface
  do i = 1, 5
    call counter(count)
  enddo
  allocate(edgeAdv%putmap(max_neigh_edges,nelemd))
  allocate(edgeAdv%getmap(max_neigh_edges,nelemd))
  allocate(edgeAdv%reverse(max_neigh_edges,nelemd))
  allocate(edgeAdv%buf(2000))
  allocate(edgeAdv%receive(2000))
  do i = 1, 5
    call counter(count)
  enddo
  !print *, tol_limiter*5D-3*5.82989072416270755D-293
  do ie = nets, nete
    do i = 1, max_neigh_edges
      edgeAdv%putmap(i, ie) = i;
      edgeAdv%getmap(i, ie) = i;
      edgeAdv%reverse(i, ie) = .true.
    enddo
  enddo
  !call biharmonic_wk_scalar(elem, Qtens_biharmonic, deriv, edgeAdv, hybrid, nets, nete)
  call edgeVpack(elem, edgeAdv, nets, nete, np1_qdp, DSSopt)

  !param_s%qdp = loc(elem(nets)%state%Qdp)
  !param_s%rkstage = rkstage
  !param_s%n0_qdp = n0_qdp
  !param_s%np1_qdp = np1_qdp
  !param_s%nets = nets
  !param_s%nete = nete
  !param_s%qsize = qsize
  !param_s%qsize_d = qsize_d
  !param_s%step_elem = (loc(elem(nets+1)%state%Qdp) - loc(elem(nets)%state%Qdp))/8
  !call athread_spawn(slave_qdp_time_avg, param_s)
  !call athread_join();
  !param_s%qdp_s_ptr = loc(elem(nets)%state%Qdp(:,:,:,:,:))
  !param_s%qdp_leap_ptr = loc(elem((nets+1))%state%Qdp(:,:,:,:,:))
  !param_s%dp_s_ptr = loc(elem(nets)%derived%dp(:,:,:))
  !param_s%divdp_proj_s_ptr = loc(elem(nets)%derived%divdp_proj(:,:,:))
  !param_s%Qtens_biharmonic = loc(Qtens_biharmonic(1,1,1,1,nets))
  !param_s%qmax = loc(qmax(1,1,nets))
  !param_s%qmin = loc(qmin(1,1,nets))
  !param_s%dt = dt
  !param_s%nets = nets
  !param_s%nete = nete
  !param_s%np1_qdp = np1_qdp
  !param_s%n0_qdp = n0_qdp
  !param_s%DSSopt = DSSopt
  !param_s%rhs_multiplier = rhs_multiplier
  !param_s%qsize = qsize
  !param_s%qsize_d = qsize_d
  !call athread_spawn(slave_euler_step, param_s)
  !call athread_join()
  !param_2d_s%qdp_s_ptr = loc(elem(nets)%state%Qdp(:,:,:,:,:))
  !param_2d_s%qdp_leap_ptr = loc(elem((nets+1))%state%Qdp(:,:,:,:,:))
  !param_2d_s%divdp_proj = loc(elem(nets)%derived%divdp_proj(:,:,:))
  !param_2d_s%dp = loc(elem(nets)%derived%dp(:,:,:))
  !param_2d_s%vn0 = loc(elem(nets)%derived%vn0(:,:,:,:))
  !param_2d_s%Dvv = loc(deriv%Dvv)
  !param_2d_s%Dinv = loc(elem(nets)%Dinv(:,:,:,:))
  !param_2d_s%metdet = loc(elem(nets)%metdet(:,:))
  !param_2d_s%rmetdet = loc(elem(nets)%rmetdet(:,:))
  !param_2d_s%Qtens_biharmonic = loc(Qtens_biharmonic(1,1,1,1,nets))
  !param_2d_s%divdp = loc(elem(nets)%derived%divdp)
  !param_2d_s%dpdiss_biharmonic = loc(elem(nets)%derived%dpdiss_biharmonic)
  !param_2d_s%spheremp = loc(elem(nets)%spheremp)
  !param_2d_s%qmax = loc(qmax(1,1,nets))
  !param_2d_s%qmin = loc(qmin(1,1,nets))
  !param_2d_s%dt = dt
  !param_2d_s%rrearth = rrearth
  !param_2d_s%nu_p = nu_p
  !param_2d_s%nu_q = nu_q
  !param_2d_s%nets = nets
  !param_2d_s%nete = nete
  !param_2d_s%rhs_multiplier = rhs_multiplier
  !param_2d_s%qsize = qsize
  !param_2d_s%qsize_d = qsize_d
  !param_2d_s%n0_qdp = n0_qdp
  !param_2d_s%np1_qdp = np1_qdp
  !param_2d_s%limiter_option = limiter_option
  !param_2d_s%rhs_viss = rhs_viss
  !call cpu_time(time_s)
  !call athread_spawn(slave_euler_v, param_2d_s)
  !call athread_join()
  !call cpu_time(time_e)
  !print *, ">>>>>>>>>>>>>time:" , time_e - time_s
  !param_compute_s%Qtens_biharmonic = loc(Qtens_biharmonic(1,1,1,1,nets))
  !param_compute_s%spheremp = loc(elem(nets)%spheremp)
  !param_compute_s%dp0 = loc(dp0)
  !param_compute_s%dt = dt;
  !param_compute_s%nu_q = nu_q
  !param_compute_s%nets = nets
  !param_compute_s%nete = nete
  !param_compute_s%qsize = qsize
  !param_compute_s%step_elem = (loc(elem(nets+1)%spheremp) - loc(elem(nets)%spheremp))/8
  !param_compute_s%rhs_viss = rhs_viss
  !call athread_spawn(slave_compute_overlap, param_compute_s)
  !call athread_join()
!
  !param_rhs_s%Qtens_biharmonic = loc(Qtens_biharmonic(1,1,1,1,nets))
  !param_rhs_s%dpdiss_ave = loc(elem(nets)%derived%dpdiss_ave)
  !param_rhs_s%dp0 = loc(dp0)
  !param_rhs_s%nets = nets
  !param_rhs_s%nete = nete
  !param_rhs_s%qsize = qsize
  !param_rhs_s%step_elem = (loc(elem(nets+1)%derived%dpdiss_ave) - loc(elem(nets)%derived%dpdiss_ave))/8
  !call athread_spawn(slave_rhs_multiplier, param_rhs_s)
  !call athread_join()

#if 0
#define PRINT_QDP
#ifdef PRINT_QDP
  do ie=nets,nete
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            !print *, elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
            print *, Qtens_temp(i,j,k,q,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
#endif
!!#define PRINT_QTEN
!#ifdef PRINT_QTEN
!  do ie=nets,nete
!    do q=1,qsize
!      do k=1,nlev
!        do j=1,np
!          do i=1,np
!            print *, dp_star_temp(i,j,k,q,ie)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!#else
!  do ie=nets,nete
!    do q=1,qsize
!      do k=1, nlev
!        print *, qmax(k,q,ie), qmin(k,q,ie)
!      enddo
!    enddo
!  enddo
!
!  qbeg = 1
!  qend = qsize
!  kbeg = 1
!  kend = nlev
!  rhs_viss = 0
!#endif
#endif

  deallocate(edgeAdv%putmap)
  deallocate(edgeAdv%getmap)
  deallocate(edgeAdv%reverse)
  deallocate(edgeAdv%buf)
  deallocate(edgeAdv%receive)
end subroutine euler_step

subroutine biharmonic_wk_scalar(elem, qtens, deriv, edgeq, hybrid, nets, nete)
use kinds , only: real_kind
use physical_constants, only: rrearth
use element_mod, only: element_t
use control_mod, only : north, south, east, west, neast, nwest, seast, swest
use dimensions_mod, only: np, npdg, nlev, qsize, max_corner_elem, max_neigh_edges, nelemd
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use edgetype_mod, only       : EdgeBuffer_t
use control_mod, only : hypervis_scaling, hypervis_power
implicit none

type(hybrid_t),  intent(in) :: hybrid
type(element_t), intent(inout), target :: elem(:)
integer :: nets, nete
real(kind=real_kind), dimension(np, np, nlev, qsize, nets:nete) :: Qtens
type(EdgeBuffer_t),  intent(inout) :: edgeq
type(derivative_t), intent(in) :: deriv


! local
!integer :: k, kptr, i, j, ie, ic, q
!real(kind=real_kind), dimension(np,np) :: lap_p
!integer :: reverse(4, nets:nete)
!integer :: step_elem
!external :: slave_biharmonic_wk_scalar
!type param_t
!  integer*8 :: Dinv, rspheremp, spheremp, variable_hyperviscosity, tensorVisc, Qtens     \
!      , Dvv, buf
!  integer*8 :: putmap, reverse
!  real(kind=real_kind) :: rrearth, hypervis_scaling, hypervis_power
!  integer :: nets, nete, qsize, step_elem
!end type param_t
!type(param_t) :: param_s
!do ie = nets, nete
!  if (edgeq%reverse(west, ie)) then
!    reverse(1, ie) = 1
!  else
!    reverse(1, ie) = 0
!  endif
!  if (edgeq%reverse(east, ie)) then
!    reverse(2, ie) = 1
!  else
!    reverse(2, ie) = 0
!  endif
!  if (edgeq%reverse(south, ie)) then
!    reverse(3, ie) = 1
!  else
!    reverse(3, ie) = 0
!  endif
!  if (edgeq%reverse(north, ie)) then
!    reverse(4, ie) = 1
!  else
!    reverse(4, ie) = 0
!  endif
!enddo
!
!param_s%Dinv = loc(elem(nets)%Dinv)
!param_s%rspheremp = loc(elem(nets)%rspheremp)
!param_s%spheremp = loc(elem(nets)%spheremp)
!param_s%variable_hyperviscosity = loc(elem(nets)%variable_hyperviscosity)
!param_s%tensorVisc = loc(elem(nets)%tensorVisc)
!param_s%Qtens = loc(Qtens)
!param_s%Dvv = loc(deriv%Dvv)
!param_s%buf = loc(edgeq%buf)
!param_s%putmap = loc(edgeq%putmap(1,nets))
!param_s%reverse = loc(reverse)
!param_s%rrearth = rrearth
!param_s%hypervis_scaling = hypervis_scaling
!param_s%hypervis_power = hypervis_power
!param_s%nets = nets
!param_s%nete = nete
!param_s%qsize = qsize
!param_s%step_elem = (loc(elem(nets+1)%Dinv) - loc(elem(nets)%Dinv))/8
!call athread_spawn(slave_biharmonic_wk_scalar, param_s)
!call athread_join()

integer :: k, kptr, i, j, ie, ic, q

external :: slave_biharmonic_wk_scalar_2
type param_2_t
  integer*8 :: Dinv, rspheremp, spheremp, variable_hyperviscosity, tensorVisc, Qtens     \
      , Dvv, receive
  integer*8 :: getmap
  real(kind=real_kind) :: rrearth, hypervis_scaling, hypervis_power
  integer :: nets, nete, qsize, step_elem
end type param_2_t
type(param_2_t) :: param_s2

param_s2%Dinv = loc(elem(nets)%Dinv)
param_s2%rspheremp = loc(elem(nets)%rspheremp)
param_s2%spheremp = loc(elem(nets)%spheremp)
param_s2%variable_hyperviscosity = loc(elem(nets)%variable_hyperviscosity)
param_s2%tensorVisc = loc(elem(nets)%tensorVisc)
param_s2%Qtens = loc(Qtens)
param_s2%Dvv = loc(deriv%Dvv)
param_s2%receive = loc(edgeq%receive)
param_s2%getmap = loc(edgeq%getmap(1,nets))
param_s2%rrearth = rrearth
param_s2%hypervis_scaling = hypervis_scaling
param_s2%hypervis_power = hypervis_power
param_s2%nets = nets
param_s2%nete = nete
param_s2%qsize = qsize
param_s2%step_elem = (loc(elem(nets+1)%Dinv) - loc(elem(nets)%Dinv))/8
call athread_spawn(slave_biharmonic_wk_scalar_2, param_s2)
call athread_join()

#if 0
  do ie = nets, 1
    do q = 1, qsize
      do k = 1, nlev
        do j = 1, np
          do i = 1, np
            !print *, elem(ie)%variable_hyperviscosity(i,j), elem(ie)%tensorVisc(i,j,1,1) \
            !    , elem(ie)%tensorVisc(i,j,1,2), elem(ie)%tensorVisc(i,j,1,2)       \
            !    , elem(ie)%tensorVisc(i,j,2,2)
            print *, Qtens(i,j,k,q,ie)
          enddo
        enddo
      enddo
    enddo
  enddo

#endif
end subroutine biharmonic_wk_scalar

subroutine neighbor_minmax(hybrid, edgeMinMax, nets, nete, min_neigh, max_neigh)
use kinds , only: real_kind
use physical_constants, only: rrearth
use element_mod, only: element_t
use control_mod, only : north, south, east, west, neast, nwest, seast, swest
use dimensions_mod, only: np, npdg, nlev, qsize, max_corner_elem, max_neigh_edges, nelemd
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use edgetype_mod, only       : EdgeBuffer_t
use control_mod, only : hypervis_scaling, hypervis_power
implicit none

type(hybrid_t),  intent(in) :: hybrid
type(EdgeBuffer_t),  intent(inout) :: edgeMinMax
integer :: nets, nete
real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

! local
integer :: ie,k,q
integer :: kblk,qblk,kptr

external :: slave_edgespack_es
type param_edgeSpack_t
  integer*8 :: qmax, qmin, buf
  integer*8 :: putmap
  integer :: nets, nete, qsize
end type param_edgeSpack_t
type(param_edgeSpack_t) :: param_edgeSpack_s

external :: slave_edgesunpack_es
type param_edgeSunpack_t
  integer*8 :: qmax, qmin, receive
  integer*8 :: getmap
  integer :: nets, nete, qsize
end type param_edgeSunpack_t
type(param_edgeSunpack_t) :: param_edgeSunpack_s

param_edgeSpack_s%qmin = loc(min_neigh)
param_edgeSpack_s%qmax = loc(max_neigh)
param_edgeSpack_s%buf = loc(edgeMinMax%buf)
param_edgeSpack_s%putmap = loc(edgeMinMax%putmap(1,nets))
param_edgeSpack_s%nets = nets
param_edgeSpack_s%nete = nete
param_edgeSpack_s%qsize = qsize
call athread_spawn(slave_edgespack_es, param_edgeSpack_s)
call athread_join()

!param_edgeSunpack_s%qmax = loc(max_neigh)
!param_edgeSunpack_s%qmin = loc(min_neigh)
!param_edgeSunpack_s%receive = loc(edgeMinMax%receive)
!param_edgeSunpack_s%getmap = loc(edgeMinMax%getmap(1,nets))
!param_edgeSunpack_s%nets = nets
!param_edgeSunpack_s%nete = nete
!param_edgeSunpack_s%qsize = qsize
!call athread_spawn(slave_edgesunpack_es, param_edgeSunpack_s)
!call athread_join()
!print *, "hello"
end subroutine neighbor_minmax


subroutine edgeVpack(elem, edgeAdvp1, nets, nete, np1_qdp, DSSopt)
use kinds , only: real_kind
use physical_constants, only: rrearth
use element_mod, only: element_t
use control_mod, only : north, south, east, west, neast, nwest, seast, swest
use dimensions_mod, only: np, npdg, nlev, qsize, max_corner_elem, max_neigh_edges, nelemd
use hybvcoord_mod, only:hvcoord_t
use hybrid_mod, only: hybrid_t
use derivative_mod, only: derivative_t
use edgetype_mod, only       : EdgeBuffer_t
use control_mod, only : hypervis_scaling, hypervis_power
implicit none

  type(EdgeBuffer_t),  intent(inout) :: edgeAdvp1
  type(element_t), intent(inout), target :: elem(:)
  integer :: nets, nete, np1_qdp, DSSopt
  integer :: reverse(4, nets:nete)

  ! local
  integer :: ie,k,q
  integer :: kblk,qblk,kptr

  external :: slave_edgevpack_es
  type param_edgeVpack_t
    integer*8 :: v1, v2, spheremp, buf
    integer*8 :: putmap, reverse
    integer :: nets, nete, qsize, step_elem
  end type param_edgeVpack_t
  type(param_edgeVpack_t) :: param_edgeVpack_s

  external :: slave_edgevunpack_es
  type param_edgeVunpack_t
    integer*8 :: v1, v2, rspheremp, receive
    integer*8 :: getmap
    integer :: nets, nete, qsize, step_elem
  end type param_edgeVunpack_t
  type(param_edgeVunpack_t) :: param_edgeVunpack_s
  !external :: slave_edgesunpack_es
  !type param_edgeSunpack_t
  !  integer*8 :: qmax, qmin, receive
  !  integer*8 :: getmap
  !  integer :: nets, nete, qsize
  !end type param_edgeSunpack_t
  !type(param_edgeSunpack_t) :: param_edgeSunpack_s


  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  do ie = nets, nete
    if (edgeAdvp1%reverse(west, ie)) then
      reverse(1, ie) = 1
    else
      reverse(1, ie) = 0
    endif
    if (edgeAdvp1%reverse(east, ie)) then
      reverse(2, ie) = 1
    else
      reverse(2, ie) = 0
    endif
    if (edgeAdvp1%reverse(south, ie)) then
      reverse(3, ie) = 1
    else
      reverse(3, ie) = 0
    endif
    if (edgeAdvp1%reverse(north, ie)) then
      reverse(4, ie) = 1
    else
      reverse(4, ie) = 0
    endif
  enddo

  param_edgeVpack_s%v1 = loc(elem(nets)%state%Qdp(1,1,1,1,np1_qdp))
  param_edgeVpack_s%spheremp = loc(elem(nets)%spheremp)
  param_edgeVpack_s%buf = loc(edgeAdvp1%buf)
  param_edgeVpack_s%putmap = loc(edgeAdvp1%putmap(1,nets))
  param_edgeVpack_s%reverse = loc(reverse)
  param_edgeVpack_s%nets = nets
  param_edgeVpack_s%nete = nete
  param_edgeVpack_s%qsize = qsize
  param_edgeVpack_s%step_elem = (loc(elem(nets+1)%spheremp) - loc(elem(nets)%spheremp))/8
  if ( DSSopt == DSSeta         ) param_edgeVpack_s%v2 = loc(elem(nets)%derived%eta_dot_dpdn)
  if ( DSSopt == DSSomega       ) param_edgeVpack_s%v2 = loc(elem(nets)%derived%omega_p)
  if ( DSSopt == DSSdiv_vdp_ave ) param_edgeVpack_s%v2 = loc(elem(nets)%derived%divdp_proj)
  call athread_spawn(slave_edgevpack_es, param_edgeVpack_s)
  call athread_join()
  param_edgeVunpack_s%v1 = loc(elem(nets)%state%Qdp(1,1,1,1,np1_qdp))
  if ( DSSopt == DSSeta         ) param_edgeVunpack_s%v2 = loc(elem(nets)%derived%eta_dot_dpdn)
  if ( DSSopt == DSSomega       ) param_edgeVunpack_s%v2 = loc(elem(nets)%derived%omega_p)
  if ( DSSopt == DSSdiv_vdp_ave ) param_edgeVunpack_s%v2 = loc(elem(nets)%derived%divdp_proj)
  param_edgeVunpack_s%rspheremp = loc(elem(nets)%rspheremp)
  param_edgeVunpack_s%receive = loc(edgeAdvp1%receive)
  param_edgeVunpack_s%getmap = loc(edgeAdvp1%getmap(1,nets))
  param_edgeVunpack_s%nets = nets
  param_edgeVunpack_s%nete = nete
  param_edgeVunpack_s%qsize = qsize
  param_edgeVunpack_s%step_elem = (loc(elem(nets+1)%spheremp) - loc(elem(nets)%spheremp))/8
  call athread_spawn(slave_edgevunpack_es, param_edgeVunpack_s)
  call athread_join()

  !param_edgeSunpack_s%qmax = loc(max_neigh)
  !param_edgeSunpack_s%qmin = loc(min_neigh)
  !param_edgeSunpack_s%receive = loc(edgeMinMax%receive)
  !param_edgeSunpack_s%getmap = loc(edgeMinMax%getmap(1,nets))
  !param_edgeSunpack_s%nets = nets
  !param_edgeSunpack_s%nete = nete
  !param_edgeSunpack_s%qsize = qsize
  !call athread_spawn(slave_edgesunpack_es, param_edgeSunpack_s)
  !call athread_join()
  !print *, "hello"
end subroutine edgeVpack
