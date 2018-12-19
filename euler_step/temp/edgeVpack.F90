
  subroutine edgeVpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)                :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    integer,              intent(in)   :: ielem
!    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,ir,ll,iptr

    integer :: is,ie,in,iw,edgeptr

    call t_startf('edgeVpack')
    is = edge%putmap(south,ielem)
    ie = edge%putmap(east,ielem)
    in = edge%putmap(north,ielem)
    iw = edge%putmap(west,ielem)
    if (edge%nlyr < (kptr+vlyr) ) then
       print *,'edge%nlyr = ',edge%nlyr
       print *,'kptr+vlyr = ',kptr+vlyr
       call haltmp('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif
!JMD    call t_adj_detailf(+2)
!dir$ ivdep
    do k=1,vlyr
       iptr = np*(kptr+k-1)
!DIR$ UNROLL(NP)
       do i=1,np
          edge%buf(iptr+ie+i)   = v(np ,i ,k) ! East
          edge%buf(iptr+is+i)   = v(i  ,1 ,k) ! South
          edge%buf(iptr+in+i)   = v(i  ,np,k) ! North
          edge%buf(iptr+iw+i)   = v(1  ,i ,k) ! West
       enddo
    enddo

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(edge%reverse(south,ielem)) then
       do k=1,vlyr
          iptr = np*(kptr+k-1)+is
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(edge%reverse(east,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+ie
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(np,i,k)
          enddo
       enddo
    endif

    if(edge%reverse(north,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+in
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(i,np,k)
          enddo
       enddo
    endif

    if(edge%reverse(west,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+iw
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(1,i,k)
          enddo
       enddo
    endif

! SWEST
    do ll=swest,swest+max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                if (iptr > size(edge%buf)) then
                   write(6, *) 'ERROR SW: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                   call abortmp('pointer bounds ERROR SW')
                end if
                edge%buf(iptr) = v(1, 1, k)
            end do
        end if
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               if (iptr > size(edge%buf)) then
                  write(6, *) 'ERROR SE: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                  call abortmp('pointer bounds ERROR SE')
               end if
               edge%buf(iptr)=v(np, 1, k)
            end do
        end if
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               if (iptr > size(edge%buf)) then
                  write(6, *) 'ERROR NE: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                  call abortmp('pointer bounds ERROR NE')
               end if
               edge%buf(iptr) = v(np, np, k)
            end do
        end if
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               if (iptr > size(edge%buf)) then
                  write(6, *) 'ERROR NW: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                  call abortmp('pointer bounds ERROR NW')
               end if
               edge%buf(iptr) = v(1, np, k)
            end do
        end if
    end do

!JMD    call t_adj_detailf(-2)

    call t_stopf('edgeVpack')
  end subroutine edgeVpack

  
  subroutine edgeVunpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (EdgeBuffer_t),         intent(in)  :: edge

    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(np,np,vlyr)
!$dir assume_aligned v:64
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem
    !type (EdgeDescriptor_t)            :: desc

    ! Local
    integer :: i,k,ll,iptr
    integer :: is,ie,in,iw
    integer :: ks,ke,kblock
    logical :: done

!JMD    call t_adj_detailf(+2)
    call t_startf('edgeVunpack')

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)

 ! ks = 1
 ! kblock = 52
 ! ke = min(vlyr,kblock)
 ! do while ((ke<=vlyr) .and. (.not. ks>vlyr))
!    if((ielem==1) .and. (iam ==1)) then
!       print *,'edgeVunpack: ks,ke ',ks,ke
!    endif

!    nb=ceiling(real(vlyr,kind=real_kind)/real(kblock,kind=real_kind))


!    if((ielem==1) .and. (iam == 1)) then
!       print *,'edgeVunpack: vlyr:= ',vlyr
!    endif

!dir$ ivdep
    do k=1,vlyr
       iptr=np*(kptr+k-1)
!DIR$ UNROLL(NP)
       do i=1,np
          v(np ,i  ,k) = v(np ,i  ,k)+edge%receive(iptr+i+ie) ! East
          v(i  ,1  ,k) = v(i  ,1  ,k)+edge%receive(iptr+i+is) ! South
          v(i  ,np ,k) = v(i  ,np ,k)+edge%receive(iptr+i+in) ! North
          v(1  ,i  ,k) = v(1  ,i  ,k)+edge%receive(iptr+i+iw) ! West
       enddo
    end do

! SWEST
    do ll=swest,swest+max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            do k=1,vlyr
                v(1  ,1 ,k)=v(1 ,1 ,k)+edge%receive((kptr+k-1)+edge%getmap(ll,ielem)+1)
            enddo
        endif
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            do k=1,vlyr
                v(np ,1 ,k)=v(np,1 ,k)+edge%receive((kptr+k-1)+edge%getmap(ll,ielem)+1)
            enddo
        endif
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            do k=1,vlyr
                v(np ,np,k)=v(np,np,k)+edge%receive((kptr+k-1)+edge%getmap(ll,ielem)+1)
            enddo
        endif
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            do k=1,vlyr
                v(1  ,np,k)=v(1 ,np,k)+edge%receive((kptr+k-1)+edge%getmap(ll,ielem)+1)
            enddo
        endif
    end do
!    ks = ke+1
!    ke = ke+kblock
!    if(ke>vlyr) then
!       ke=vlyr
!    endif
! enddo
    call t_stopf('edgeVunpack')

!JMD    call t_adj_detailf(-2)

  end subroutine edgeVunpack
