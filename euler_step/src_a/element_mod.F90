module element_mod
use kinds, only: real_kind
use dimensions_mod, only: np, nlev, nlevp, qsize_d, qsize
  implicit none
  integer, public, parameter :: timelevels = 3

  type, public :: elem_state_t
    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)            ! velocity                           1
    real (kind=real_kind) :: T   (np,np,nlev,timelevels)              ! temperature                        2
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)              ! delta p on levels                  8
    real (kind=real_kind) :: lnps(np,np,timelevels)                   ! log surface pressure               3
    real (kind=real_kind) :: ps_v(np,np,timelevels)                   ! surface pressure                   4
    real (kind=real_kind) :: phis(np,np)                              ! surface geopotential (prescribed)  5
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                 ! Tracer concentration               6
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)               ! Tracer mass                        7
  end type elem_state_t

  type, public :: derived_state_t

    ! diagnostic variables for preqx solver

    ! storage for subcycling tracers/dynamics
    ! if (compute_mean_flux==1) vn0=time_avg(U*dp) else vn0=U at tracer-time t

    real (kind=real_kind) :: vn0  (np,np,2,nlev)                      ! velocity for SE tracer advection
    real (kind=real_kind) :: vstar(np,np,2,nlev)                      ! velocity on Lagrangian surfaces
    real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)            ! mean dp dissipation tendency, if nu_p>0
    real (kind=real_kind) :: dpdiss_ave(np,np,nlev)                   ! mean dp used to compute psdiss_tens

    ! diagnostics for explicit timestep
    real (kind=real_kind) :: phi(np,np,nlev)                          ! geopotential
    real (kind=real_kind) :: omega_p(np,np,nlev)                      ! vertical tendency (derived)
    real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)                ! mean vertical flux from dynamics

    ! semi-implicit diagnostics: computed in explict-component, reused in Helmholtz-component.
    real (kind=real_kind) :: grad_lnps(np,np,2)                       ! gradient of log surface pressure
    real (kind=real_kind) :: zeta(np,np,nlev)                         ! relative vorticity
    real (kind=real_kind) :: div(np,np,nlev,timelevels)               ! divergence

    ! tracer advection fields used for consistency and limiters
    real (kind=real_kind) :: dp(np,np,nlev)                           ! for dp_tracers at physics timestep
    real (kind=real_kind) :: divdp(np,np,nlev)                        ! divergence of dp
    real (kind=real_kind) :: divdp_proj(np,np,nlev)                   ! DSSed divdp

    ! forcing terms for CAM
    real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, 1)                ! tracer forcing
    real (kind=real_kind) :: FM(np,np,2,nlev, 1)                      ! momentum forcing
    real (kind=real_kind) :: FT(np,np,nlev, 1)                        ! temperature forcing
    real (kind=real_kind) :: omega_prescribed(np,np,nlev)             ! prescribed vertical tendency

    ! forcing terms for both CAM and HOMME
    ! FQps for conserving dry mass in the presence of precipitation

    real (kind=real_kind) :: pecnd(np,np,nlev)                        ! pressure perturbation from condensate
    real (kind=real_kind) :: FQps(np,np,timelevels)                   ! forcing of FQ on ps_v
  end type derived_state_t

  type, public :: element_t
     type (elem_state_t) :: state
     type (derived_state_t)   :: derived
     real (kind=real_kind)    :: met(np,np,2,2)                       ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metinv(np,np,2,2)                    ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metdet(np,np)                        ! g = SQRT(det(g_ij)) on velocity and pressure grid
     real (kind=real_kind)    :: rmetdet(np,np)                       ! 1/metdet on velocity pressure grid
     real (kind=real_kind)    :: D(np,np,2,2)                         ! Map covariant field on cube to vector field on the sphere
     real (kind=real_kind)    :: Dinv(np,np,2,2)                      ! Map vector field
     real (kind=real_kind)    :: spheremp(np,np)                      ! mass matrix on v and p grid
     real (kind=real_kind)    :: rspheremp(np,np)                     ! inverse mass matrix on v and p grid
     real (kind=real_kind)    :: variable_hyperviscosity(np,np)       ! hyperviscosity based on above
     real (kind=real_kind)    :: hv_courant                           ! hyperviscosity courant number
     real (kind=real_kind)    :: tensorVisc(np,np,2,2)                !og, matrix V for tensor viscosity

  end type element_t

end module element_mod
