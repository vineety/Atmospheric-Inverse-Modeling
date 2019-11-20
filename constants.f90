    module constants
    implicit none
    integer, parameter:: rkind=kind(1.d0) !selected_real_kind(15,307) ! precision of real numbers can only be 4 or 8
    !integer, parameter:: rkind=4 !for single precision
    integer, parameter:: ikind=4 ! precision of integers can only be 4 or 8
    character(len=6),parameter::accuracy='double'! single for single precision , double for double precision
    !character(len=6),parameter::accuracy='single'! single for single precision , double for double precision
    character(len=9),parameter:: mpi_dtype='mpi_real8'
    !integer(kind=ikind):: myid=1
    !integer(kind=ikind):: numprocs=1
    end module constants
