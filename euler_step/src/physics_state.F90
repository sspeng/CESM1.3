module physics_types
!-------------------------------------------------------------------------------
  type physics_state
     integer                                     :: &
          lchnk,                &! chunk index
          ngrdcol,              &! -- Grid        -- number of active columns (on the grid)
          psetcols=0,           &! --             -- max number of columns set - if subcols = pcols*psubcols, else = pcols
          ncol=0                 ! --             -- sum of nsubcol for all ngrdcols - number of active columns
     real(r8), dimension(:), allocatable         :: &
          lat,     &! latitude (radians)
          lon,     &! longitude (radians)
          ps,      &! surface pressure
          psdry,   &! dry surface pressure
          phis,    &! surface geopotential
          ulat,    &! unique latitudes  (radians)
          ulon      ! unique longitudes (radians)
     real(r8), dimension(:,:),allocatable        :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s)
          pmid,    &! midpoint pressure (Pa)
          pmiddry, &! midpoint pressure dry (Pa)
          pdel,    &! layer thickness (Pa)
          pdeldry, &! layer thickness dry (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          rpdeldry,&! recipricol layer thickness dry (Pa)
          lnpmid,  &! ln(pmid)
          lnpmiddry,&! log midpoint pressure dry (Pa)
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
          zm        ! geopotential height above surface at midpoints (m)

     real(r8), dimension(:,:,:),allocatable      :: &
          q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

     real(r8), dimension(:,:),allocatable        :: &
          pint,    &! interface pressure (Pa)
          pintdry, &! interface pressure dry (Pa)
          lnpint,  &! ln(pint)
          lnpintdry,&! log interface pressure dry (Pa)
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(:),allocatable          :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     integer, dimension(:),allocatable           :: &
          latmapback, &! map from column to unique lat for that column
          lonmapback, &! map from column to unique lon for that column
          cid        ! unique column id
     integer :: ulatcnt, &! number of unique lats in chunk
                uloncnt   ! number of unique lons in chunk

  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend

     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

     real(r8), dimension(:,:),allocatable        :: dtdt, dudt, dvdt
     real(r8), dimension(:),  allocatable        :: flx_net
     real(r8), dimension(:),  allocatable        :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend

     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls = .false.,               &! true if dsdt is returned
          lu = .false.,               &! true if dudt is returned
          lv = .false.                 ! true if dvdt is returned

     logical,dimension(pcnst) ::  lq = .false.  ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(:,:),allocatable   :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(:,:,:),allocatable :: &
          q                  ! consituent tendencies (kg/kg/s)

! boundary fluxes
     real(r8), dimension(:),allocatable     ::&
          hflux_srf,     &! net heat flux at surface (W/m2)
          hflux_top,     &! net heat flux at top of model (W/m2)
          taux_srf,      &! net zonal stress at surface (Pa)
          taux_top,      &! net zonal stress at top of model (Pa)
          tauy_srf,      &! net meridional stress at surface (Pa)
          tauy_top        ! net meridional stress at top of model (Pa)
     real(r8), dimension(:,:),allocatable   ::&
          cflx_srf,      &! constituent flux at surface (kg/m2/s)
          cflx_top        ! constituent flux top of model (kg/m2/s)

  end type physics_ptend

  subroutine cnst_get_ind (name, ind, abort)
!-----------------------------------------------------------------------
!
! Purpose: Get the index of a constituent
!
! Author:  B.A. Boville
!
!-----------------------------Arguments---------------------------------
!
    character(len=*),  intent(in)  :: name  ! constituent name
    integer,           intent(out) :: ind   ! global constituent index (in q array)
    logical, optional, intent(in)  :: abort ! optional flag controlling abort

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    logical :: abort_on_error
!-----------------------------------------------------------------------

! Find tracer name in list
    do m = 1, pcnst
       if (name == cnst_name(m)) then
          ind  = m
          return
       end if
    end do

! Unrecognized name
    abort_on_error = .true.
    if ( present(abort) ) abort_on_error = abort

    if ( abort_on_error ) then
       write(iulog,*) 'CNST_GET_IND, name:', name,  ' not found in list:', cnst_name(:)
       call endrun('CNST_GET_IND: name not found')
    end if

! error return
    ind = -1

  end subroutine cnst_get_ind

    subroutine state_cnst_min_nz(lim, qix, numix)
      ! Small utility function for setting minimum nonzero
      ! constituent concentrations.

      ! Lower limit and constituent index
      real(r8), intent(in) :: lim
      integer,  intent(in) :: qix
      ! Number concentration that goes with qix.
      ! Ignored if <= 0 (and therefore constituent is not present).
      integer,  intent(in) :: numix

      if (numix > 0) then
         ! Where q is too small, zero mass and number
         ! concentration.
         where (state%q(:ncol,:,qix) < lim)
            state%q(:ncol,:,qix) = 0._r8
            state%q(:ncol,:,numix) = 0._r8
         end where
      else
         ! If no number index, just do mass.
          where (state%q(:ncol,:,qix) < lim)
             state%q(:ncol,:,qix) = 0._r8
          end where
      end if

    end subroutine state_cnst_min_nz
end module physics_types
