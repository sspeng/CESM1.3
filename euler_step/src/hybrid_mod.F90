module hybrid_mod
implicit none
  integer, parameter :: ncomponents = 1

  type, public :: parallel_t
    integer :: rank                       ! local rank
    integer :: root                       ! local root
    integer :: nprocs                     ! number of processes in group
    integer :: comm                       ! local communicator
    integer :: intercomm(0:ncomponents-1) ! inter communicator list
    integer :: intracomm                  ! intra node communicator
    integer :: intracommsize              ! number of MPI ranks in intra communicator
    integer :: intracommrank              ! rank in intra communicator
    integer :: commGraphFull
    integer :: commGraphInter
    integer :: commGraphIntra
    integer :: groupGraphFull
    logical :: masterproc
  end type

  type, public :: hybrid_t
    type (parallel_t) :: par
    integer           :: ithr
    integer           :: localsense
    integer           :: nthreads
    integer           :: ibeg, iend
    integer           :: kbeg, kend
    integer           :: qbeg, qend
    logical           :: masterthread
  end type

end module hybrid_mod
