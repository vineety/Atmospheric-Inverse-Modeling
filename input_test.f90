    program input_test
    !***********************************************************************
    !
    !written by : Mae Qiu
    !date : 7/10/2013
    !purpose : test whether the input files for the inverse modeling code
    !          are in the correct format
    !***********typset and definitions**************************************
    use file_io
    use constants
    use initialize
    use library_generic
    implicit none
    real(kind=rkind),allocatable,dimension(:,:)::h
    real(kind=rkind),allocatable,dimension(:,:)::x
    real(kind=rkind),allocatable,dimension(:,:)::sp_coord
    real(kind=rkind),allocatable,dimension(:,:)::parameters
    real(kind=rkind),allocatable,dimension(:,:)::r
    real(kind=rkind),allocatable,dimension(:)::z
    real(kind=rkind),allocatable,dimension(:)::sp
    real(kind=rkind),allocatable,dimension(:)::q_diag
    integer(kind=ikind)::nargin
    integer(kind=ikind)::i
    integer(kind=ikind)::nnz_h
    integer(kind=ikind)::nnz_x
    character(len=200)::argin
    character(len=200)::fname
    character(len=200)::config_file ! character to hold file name
    logical::dir_e
    type(param)::prm
    !**********execution part**************************************************
    nargin=iargc()
    if (nargin .ne. 1) then
       write(*,'(a)') 'usage: ./input_test.exe config_file'
       write(*,'(a)') 'enter the name of the config file:'
       read(*,*) config_file
    else
       call getarg(1,config_file)
    endif
    config_file=adjustl(trim(config_file))
    prm=set_all_args(config_file)

    fname=trim(adjustl(prm%inp_direc))
    inquire(directory=fname, exist=dir_e)
    if ( .not. dir_e ) then
       write(*,*) 'Input directory=',fname
       write(*,*) "Input directory doesn't exist!"
       stop
    endif

    fname=trim(adjustl(prm%out_direc))
    inquire(directory=fname, exist=dir_e )
    if ( .not. dir_e ) then
       write(*,*) 'Output directory=',fname
       write(*,*) "Output directory doesn't exist!"
       stop
    endif

    if (prm%invtype /='geostat' .and. prm%invtype /='bayesian') then
       write(*,*) 'Inversion type can only be geostat or bayesian!'
       stop
    endif
       
    if (prm%covtype /='space_time' .and. prm%covtype /='spatial' &
         .and. prm%covtype /='independent') then
       write(*,*) 'Covariance type can only be space_time, spatial or independent!'
       stop
    endif
       
    if (prm%hincore /='yes' .and. prm%hincore /='no') then
       write(*,*) 'H_incore can only be yes or no!'
       stop
    endif
       
    if (prm%xincore /='yes' .and. prm%xincore /='no') then
       write(*,*) 'X_incore can only be yes or no!'
       stop
    endif
       
    if (prm%xsparse /='yes' .and. prm%xsparse /='no') then
       write(*,*) 'X_sparse can only be yes or no!'
       stop
    endif
    
    if (prm%locations<=0) then
       write(*,*) 'Number_flux_locations have to be greater than 0!'
       stop
    endif

    if (prm%observations<=0) then
       write(*,*) 'Number_observations have to be greater than 0!'
       stop
    endif

    if (mod(24,prm%withinday) /=0) then
       write(*,*) 'Check the Time_resolution.'
       write(*,*) 'For hourly data, it can only be 1,2,3,4,6,8,12 hourly'
       stop
    endif

    if (prm%overlap <0) then
       write(*,*) 'Days_discard has to be greater than or equal to 0!'
       stop
    endif

    write(*,*) '***********************************************************'
800 format(a,f10.3,a,i0,a,i0,a,f10.3)
900 format(a,f10.3,a,f10.3)
    call read_binary_matrix(prm%inp_direc,prm%filename_lat_long,sp_coord)
    write(*,"(a,2f10.3)") '1st pair in spatial coordinates: ',sp_coord(1,:)
    write(*,"(a,2f10.3)") 'last pair in spatial coordinates: ',sp_coord(size(sp_coord,1),:)

    if (prm%covtype =='independent') then
       call read_binary_vector(prm%inp_direc,prm%filename_parameters,q_diag)
       write(*,900) '1st element in Q=',q_diag(1),', last element in Q=',q_diag(size(q_diag,1))
    else
       call read_binary_matrix(prm%inp_direc,prm%filename_parameters,parameters)
       prm%variance=parameters(:,1)
       prm%sclength=parameters(:,2)
       prm%tclengthd=parameters(:,3)
       write(*,900) '1st variance=',prm%variance(1),', last variance=',&
                prm%variance(prm%size)
       write(*,900) '1st sclength=',prm%sclength(1),', last sclength=',&
                prm%sclength(prm%size)
       write(*,900) '1st tclengthd=',prm%tclengthd(1),', last tclengthd=',&
                prm%tclengthd(prm%size)
    endif
    write(*,*) '--------------------------------------------------------'
    write(*,*) '-----------------Checking H & x strips:-----------------'

    allocate(h(prm%observations,prm%locations))
    if (prm%invtype=='geostat') then
       allocate(x(prm%locations,prm%covariates))
    endif
    do i=1,sum(prm%days)*prm%withinday
       call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,i,&
                               prm%observations,prm%locations,h,nnz_h)
       write(*,"(a,i0,a,i0)") 'H strip #',i,', nnz=',nnz_h
       if (nnz_h .ne. 0) then
          write(*,800) 'h(1,1)=',h(1,1),',   h(',prm%observations,&
               ',',prm%locations,')=',h(prm%observations,prm%locations)
       endif
       if (prm%invtype=='geostat') then
          call sparse_full_return(prm%inp_direc,prm%name_str_x,prm%ext,i,&
                               prm%locations,prm%covariates,x,nnz_x)
          write(*,"(a,i0,a,i0)") 'x strip #',i,', nnz=',nnz_x
          if (nnz_x .ne. 0) then
             write(*,800) 'x(1,1)=',x(1,1),',   x(',prm%locations,&
               ',',prm%covariates,')=',h(prm%locations,prm%covariates)
          endif
       endif
       write(*,*) '--------------------------------------------------------'
    enddo

    write(*,*) 'Checking r matrix:'
    fname=prm%filename_r
    call read_binary_vector(prm%inp_direc,fname,r)
    write(*,900) '1st element in r=',r(1),', last element in r=',r(size(r,1))
    !call read_binary_matrix(prm%inp_direc,fname,r)
    !write(*,900) '1st element in r=',r(1,1),', last element in r=',r(size(r,1),size(r,1))
    write(*,*) '--------------------------------------------------------'

    write(*,*) 'Checking z vector:'
    fname=prm%filename_z
    call read_binary_vector(prm%inp_direc,fname,z)
    write(*,900) '1st element in z=',z(1),', last element in z=',z(size(z,1))
    write(*,*) '--------------------------------------------------------'

    if (prm%invtype =='bayesian') then
       fname=prm%filename_sp
       call read_binary_vector(prm%inp_direc,fname,sp)
       write(*,*) 'Checking sp vecotr:'
       write(*,900) '1st element in sp=',sp(1),', last element in sp=',sp(size(sp,1))
    endif

    write(*,*) '***********************************************************'
    write(*,*) 'Please make sure the above 1st and last elements for the'
    write(*,*) 'corresponding matrix/vector are the same as your input data.'
    write(*,*) 'Once you''re sure that the input data is read in correctly and'
    write(*,*) 'this test program doesn''t give you any other error messages, '
    write(*,*) 'you are ready to run the inverse modeling code!'
    write(*,*) '***********************************************************'
    end program input_test
