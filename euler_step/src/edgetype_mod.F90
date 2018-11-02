module edgetype_mod

  use kinds, only : int_kind, long_kind, log_kind, real_kind
  implicit none
  integer, public :: initedgebuffer_callid = 0
  type, public :: EdgeBuffer_t
     real (kind=real_kind), dimension(:), allocatable :: buf
     real (kind=real_kind), dimension(:), allocatable :: receive
     integer(kind=int_kind), pointer :: putmap(:,:) => null()
     integer(kind=int_kind), pointer :: getmap(:,:) => null()
     logical(kind=log_kind), pointer :: reverse(:,:) => null()
     integer(kind=int_kind), pointer :: moveLength(:) => null()
     integer(kind=int_kind), pointer :: movePtr(:) => null()

     integer(kind=int_kind), pointer :: rcountsFull(:) => null()
     integer(kind=int_kind), pointer :: scountsFull(:) => null()
     integer(kind=int_kind), pointer :: sdisplsFull(:) => null()
     integer(kind=int_kind), pointer :: rdisplsFull(:) => null()

     integer(kind=int_kind), pointer :: rcountsInter(:) => null()
     integer(kind=int_kind), pointer :: scountsInter(:) => null()
     integer(kind=int_kind), pointer :: sdisplsInter(:) => null()
     integer(kind=int_kind), pointer :: rdisplsInter(:) => null()

     integer(kind=int_kind), pointer :: rcountsIntra(:) => null()
     integer(kind=int_kind), pointer :: scountsIntra(:) => null()
     integer(kind=int_kind), pointer :: sdisplsIntra(:) => null()
     integer(kind=int_kind), pointer :: rdisplsIntra(:) => null()

     integer(kind=int_kind), pointer :: getDisplsFull(:) => null()
     integer(kind=int_kind), pointer :: putDisplsFull(:) => null()

     integer(kind=int_kind), dimension(:), allocatable :: Rrequest,Srequest
     integer(kind=int_kind), dimension(:,:), allocatable :: status
     integer :: nlyr ! Number of layers
     integer :: nbuf ! total size of message passing buffer, includes vertical levels
     integer :: nInter, nIntra
     integer :: id
     integer :: bndry_type
     integer :: tag
     integer :: win
     integer(kind=long_kind) :: winSize
  end type EdgeBuffer_t
end module
