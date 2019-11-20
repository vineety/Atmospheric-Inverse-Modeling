    module initialize
    !****************************************************************************
    !  modul: initialize
    !  purpose: initialize the parameters from a config file
    !****************************************************************************
    use constants
    use file_io
    implicit none

    type param
        character(len=100)::inp_direc
        character(len=100)::out_direc
        character(len=100)::filename_lat_long
        character(len=100)::filename_r
        character(len=100)::filename_z
        character(len=100)::filename_sp
        character(len=100)::filename_parameters
        character(len=15)::invtype
        character(len=15)::covtype
        character(len=15)::hincore
        character(len=15)::xincore
        character(len=15)::xsparse
        character(len=200)::name_str_h
        character(len=200)::name_str_x
        character(len=200)::ext
        integer(kind=ikind)::size
        integer(kind=ikind)::locations
        integer(kind=ikind)::observations
        integer(kind=ikind)::withinday
        integer(kind=ikind)::overlap
        integer(kind=ikind)::covariates
        integer(kind=ikind),allocatable,dimension(:)::days
        real(kind=rkind)::time_resolution
        real(kind=rkind),allocatable, dimension(:)::variance
        real(kind=rkind),allocatable, dimension(:)::sclength
        real(kind=rkind),allocatable, dimension(:)::tclengthd
        real(kind=rkind),allocatable, dimension(:)::tclengthwd
    end type param

    contains

    function set_all_args (file_name) result(args)
    ! purpose:
    ! set all the parameters read from a config file
    !%%%%%%%%%%%%%%%%%% Definition Part %%%%%%%%%%%%%%%%%%%%%%%%%
    implicit none
    character(len=*) file_name
    character(len=200) arg_type
    character(len=200) fname
    character(len=6) repeat
    integer(kind=ikind):: configerr
    integer(kind=ikind):: tmp
    integer(kind=ikind) :: nnz 
    integer(kind=ikind),allocatable, dimension(:,:)::coord
    real(kind=rkind),allocatable,dimension(:):: values
    type(param)::args
    !%%%%%%%%%%%%%%%%%%%%%%% open file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    open(unit=1,file=trim(adjustl(file_name)),form="formatted",&
    iostat=configerr,status="old",action="read")
    if (configerr /=0) then
       write(*,*) 'Error opening the config file: ',trim(adjustl(file_name))
       stop
    endif

    !%%%%%%%%%%%%%%%%%%%%%%% read inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    read(1,*) arg_type, args%inp_direc
    args%inp_direc=trim(adjustl(args%inp_direc))
    read(1,*) arg_type, args%out_direc
    args%out_direc=trim(adjustl(args%out_direc))
    read(1,*) arg_type, args%invtype
    args%invtype=trim(adjustl(args%invtype))
    read(1,*) arg_type,args%covtype
    args%covtype=trim(adjustl(args%covtype))
    read(1,*) arg_type, args%hincore
    args%hincore=trim(adjustl(args%hincore))
    read(1,*) arg_type, args%xincore
    args%xincore=trim(adjustl(args%xincore))
    read(1,*) arg_type, args%xsparse
    args%xsparse=trim(adjustl(args%xsparse))
    read(1,*) arg_type, args%filename_lat_long
    args%filename_lat_long=trim(adjustl(args%filename_lat_long))
    read(1,*) arg_type, args%filename_r
    args%filename_r=trim(adjustl(args%filename_r))
    read(1,*) arg_type, args%filename_z
    args%filename_z=trim(adjustl(args%filename_z))
    read(1,*) arg_type, args%filename_sp
    args%filename_sp=trim(adjustl(args%filename_sp))
    read(1,*) arg_type, args%filename_parameters
    args%filename_parameters=trim(adjustl(args%filename_parameters))
    read(1,*) arg_type, args%name_str_h
    args%name_str_h=trim(adjustl(args%name_str_h))
    read(1,*) arg_type, args%name_str_x
    args%name_str_x=trim(adjustl(args%name_str_x))
    read(1,*) arg_type, args%ext
    args%ext=trim(adjustl(args%ext))
    !%%%%%%%%%%% read in number of covariates from X strips %%%%%%%% 
    if (args%xsparse=='yes') then    
       fname=trim(adjustl(args%name_str_x))//'1'//trim(adjustl(args%ext))
       call sparse_read_coo(args%inp_direc,fname,nnz,coord,values)
       !the last entry of coord is the number of rows & cols of X1 matrix
       !and the number of covariates is equal to the # of cols in X1 matrix
       args%covariates=coord(nnz,2)  
    elseif (args%xsparse=='no') then
       fname=trim(adjustl(args%name_str_x))//'f'//trim(adjustl(args%ext))
       call read_rows_cols(args%inp_direc,fname,tmp,args%covariates)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    read(1,*) arg_type, args%locations
    read(1,*) arg_type, args%observations
    read(1,*) arg_type, args%size
    !%%%%%%%%%%%%%%%%%%%%%%% allocate arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(args%variance(args%size),args%sclength(args%size),&
    args%tclengthd(args%size),args%tclengthwd(args%size),args%days(args%size))
    !%%%%%%%%%%%%%%%%%%%%% read in values for arrays%%%%%%%%%%%%%%%%%%%%%
    read(1,*,iostat=configerr) arg_type, args%days
    if (configerr /=0) then
       write(*,*) 'Check the input for days.'
       write(*,*) 'The number of elements in days has to be equal to the size'
       write(*,*) 'given at the first line of the config file.'
       stop
    endif
    read(1,*) arg_type, args%time_resolution
    if (args%time_resolution>=1) then
       args%withinday=1
    else
       args%withinday=nint(1/args%time_resolution)
    endif
    read(1,*) arg_type, args%overlap
    !%%%%%%%%%%%%%%%%%%%%% write read in values %%%%%%%%%%%%%%%%%%%%%%%%
    write(repeat,'(i6)') args%size
    write(*,"(a,a)") 'Input Directory: ' ,trim(args%inp_direc)
    write(*,"(a,a)") 'Output Directory: ' ,trim(args%out_direc)
    write(*,"(a,a)") 'Filename_lat_long: ' ,trim(args%filename_lat_long)
    write(*,"(a,a)") 'R Vector File Name: ' ,trim(args%filename_r)
    write(*,"(a,a)") 'Z Vector File Name: ' ,trim(args%filename_z)
    write(*,"(a,a)") 'S Prior Vector File Name: ' ,trim(args%filename_sp)
    write(*,"(a,a)") 'Q Parameter/Diagonal File Name: ' ,trim(args%filename_parameters)
    write(*,"(a,a)") 'Inversion type: ' ,trim(args%invtype)
    write(*,"(a,a)") 'Covariance type: ' ,trim(args%covtype)
    write(*,"(a,a)") 'H_incore: ' ,trim(args%hincore)
    write(*,"(a,a)") 'X_incore: ' ,trim(args%xincore)
    write(*,"(a,a)") 'X_sparse: ' ,trim(args%xsparse)
    write(*,"(a,i0)") 'Number of flux locations: ' ,args%locations
    write(*,"(a,i0)") 'Number of Observations: ', args%observations
    write(*,"(a,i0)") 'Total Covariates: ', args%covariates
    write(*,"(a,i0)") 'Time arregation periods: ' ,args%size
    write(*,"(a,"//trim(adjustl(repeat))//"I4)") 'Time arregation: ', args%days
    write(*,"(a,f6.3)") 'Time resolution: ', args%time_resolution
    write(*,"(a,i0)") 'Days discard: ', args%overlap
    !%%%%%%%%%%%%%%%%%% add Days discard to the 1st time period %%%%%%%%%%%%%%%%
    args%days(1)=args%days(1)+args%overlap
    close(1)
    end function set_all_args
    end module initialize
