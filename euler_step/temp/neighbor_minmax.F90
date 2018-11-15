#ifdef _PRIM
subroutine neighbor_minmax(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
!
! compute Q min&max over the element and all its neighbors
!
!
type (hybrid_t)      , intent(in) :: hybrid
type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
integer :: nets, nete
real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

! local
integer :: ie,k,q
integer :: kblk,qblk,kptr
!real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
!real (kind=real_kind) :: Qmax(np,np,nlev,qsize)


   kblk = nlev  ! calculate size of the block of vertical levels
   qblk = qsize   ! calculate size of the block of tracers

   do ie=nets,nete
      do q = 1, qsize
         kptr = nlev*(q - 1)
         call  edgeSpack(edgeMinMax,min_neigh(:,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1)
         call  edgeSpack(edgeMinMax,max_neigh(:,q,ie),kblk,kptr,ie)
      enddo
   enddo

   call bndry_exchangeS(hybrid,edgeMinMax,location='neighbor_minmax')

   do ie=nets,nete
      do q=1,qsize
         kptr = nlev*(q - 1)
         call  edgeSunpackMIN(edgeMinMax,min_neigh(:,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1)
         call  edgeSunpackMAX(edgeMinMax,max_neigh(:,q,ie),kblk,kptr,ie)
         do k=1,nlev
            min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
         enddo
      enddo
   enddo

end subroutine neighbor_minmax

subroutine edgeSunpackMAX(edge,v,vlyr,kptr,ielem)
  use dimensions_mod, only : np, max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest

  type (EdgeBuffer_t),         intent(in)  :: edge
  integer,               intent(in)  :: vlyr
  real (kind=real_kind), intent(inout) :: v(vlyr)
  integer,               intent(in)  :: kptr
  integer,               intent(in)  :: ielem


  ! Local

  integer :: i,k,l,iptr
  integer :: is,ie,in,iw

  threadsafe=.false.

  is=edge%getmap(south,ielem)
  ie=edge%getmap(east,ielem)
  in=edge%getmap(north,ielem)
  iw=edge%getmap(west,ielem)
  do k=1,vlyr
     iptr=(kptr+k-1)
     v(k) = MAX(v(k),edge%receive(iptr+is+1),edge%receive(iptr+ie+1),edge%receive(iptr+in+1),edge%receive(iptr+iw+1))
  end do

! SWEST
  do l=swest,swest+max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MAX(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MAX(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MAX(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MAX(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

end subroutine edgeSunpackMAX

subroutine edgeSunpackMIN(edge,v,vlyr,kptr,ielem)
  use dimensions_mod, only : np, max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest

  type (EdgeBuffer_t),         intent(in)  :: edge
  integer,               intent(in)  :: vlyr
  real (kind=real_kind), intent(inout) :: v(vlyr)
  integer,               intent(in)  :: kptr
  integer,               intent(in)  :: ielem


  ! Local

  integer :: i,k,l,iptr
  integer :: is,ie,in,iw

  threadsafe=.false.

  is=edge%getmap(south,ielem)
  ie=edge%getmap(east,ielem)
  in=edge%getmap(north,ielem)
  iw=edge%getmap(west,ielem)
  do k=1,vlyr
     iptr=(kptr+k-1)
     v(k) = MIN(v(k),edge%receive(iptr+is+1),edge%receive(iptr+ie+1),edge%receive(iptr+in+1),edge%receive(iptr+iw+1))
  end do

! SWEST
  do l=swest,swest+max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MiN(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MIN(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MIN(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(edge%getmap(l,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              v(k)=MIN(v(k),edge%receive(iptr+edge%getmap(l,ielem)+1))
          enddo
      endif
  end do

end subroutine edgeSunpackMIN



subroutine edgeSpack_r8(edge,v,vlyr,kptr,ielem)
  use dimensions_mod, only : np, max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest

  type (EdgeBuffer_t)                :: edge
  integer,              intent(in)   :: vlyr
  real (kind=real_kind),intent(in)   :: v(vlyr)
  integer,              intent(in)   :: kptr
  integer,              intent(in)   :: ielem
!    type (EdgeDescriptor_t),intent(in) :: desc

  ! Local variables
  integer :: i,k,ir,ll,iptr

  integer :: is,ie,in,iw
  real (kind=real_kind) :: tmp

!    call t_adj_detailf(+2)

  is = edge%putmap(south,ielem)
  ie = edge%putmap(east,ielem)
  in = edge%putmap(north,ielem)
  iw = edge%putmap(west,ielem)
  if (edge%nlyr < (kptr+vlyr) ) then
     call haltmp('edgeSpack: Buffer overflow: size of the vertical dimension must be increased!')
  endif

  do k=1,vlyr
     iptr = kptr+k-1
     edge%buf(iptr+ie+1)   = v(k) ! East
     edge%buf(iptr+is+1)   = v(k) ! South
     edge%buf(iptr+in+1)   = v(k) ! North
     edge%buf(iptr+iw+1)   = v(k) ! West
  enddo

! SWEST
  do ll=swest,swest+max_corner_elem-1
      if (edge%putmap(ll,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              edge%buf(iptr+edge%putmap(ll,ielem)+1)=v(k)
          end do
      end if
  end do

! SEAST
  do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if (edge%putmap(ll,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              edge%buf(iptr+edge%putmap(ll,ielem)+1)=v(k)
          end do
      end if
  end do

! NEAST
  do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if (edge%putmap(ll,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              edge%buf(iptr+edge%putmap(ll,ielem)+1)=v(k)
          end do
      end if
  end do

! NWEST
  do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if (edge%putmap(ll,ielem) /= -1) then
          do k=1,vlyr
              iptr = (kptr+k-1)
              edge%buf(iptr+edge%putmap(ll,ielem)+1)=v(k)
          end do
      end if
  end do

!    call t_adj_detailf(-2)

end subroutine edgeSpack_r8
