    program inverse_modeling_main
    !***********************************************************************
    !
    !written by : vineet yadav
    !date : 7/09/2012
    !purpose : inverse modeling

    !***********typset and definitions**************************************
    !use omp_lib
    use file_io
    use constants
    use library_generic
    use library_inverse
    use parallel
    use initialize
    use steps_mod
    use distance
    implicit none
    include 'mpif.h'
    real(kind=rkind),allocatable,dimension(:,:)::hqht
    real(kind=rkind),allocatable,dimension(:,:)::sens_mat
    real(kind=rkind),allocatable,dimension(:,:)::inv_bmi
    real(kind=rkind),allocatable,dimension(:,:)::hx
    real(kind=rkind),allocatable,dimension(:,:)::sp_coord
    real(kind=rkind),allocatable,dimension(:,:)::hEff
    real(kind=rkind),allocatable,dimension(:,:)::parameters
    real(kind=rkind),allocatable,dimension(:,:)::temporal
    real(kind=rkind),allocatable,dimension(:,:)::q
    real(kind=rkind),allocatable,dimension(:,:)::junk
    real(kind=rkind),allocatable,dimension(:)::q_diag
    real(kind=rkind),allocatable,dimension(:)::Q_sum
    real(kind=rkind),allocatable,dimension(:)::z
    real(kind=rkind),allocatable,dimension(:)::r
    real(kind=rkind),allocatable,dimension(:)::sp
    real(kind=rkind),allocatable,dimension(:)::hsp
    real(kind=rkind),allocatable,dimension(:)::qhts
    real(kind=rkind), allocatable, dimension(:,:):: x
    real(kind=rkind), allocatable, dimension(:,:):: y
    real(kind=rkind), allocatable, dimension(:,:):: h
    real(kind=rkind), allocatable, dimension(:):: a
    real(kind=rkind),allocatable,dimension(:,:)::inv_psi
    real(kind=rkind), allocatable, dimension(:) :: shat1, shat2, shat
    real(kind=rkind), allocatable, dimension (:,:):: lambda_s, hqs, xs , ms
    real(kind=rkind), allocatable, dimension(:,:) :: temporary
    real(kind=rkind), allocatable, dimension(:,:) :: temporary2
    !dec$ attributes align: 32 :: hqht,hx,sp_coord
    !dec$ attributes align: 32 :: inv_bmi,hEff,temporal
    !dec$ attributes align: 32 :: shat,Q_sum,z,r,parameters
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind)::stime, etime
    integer(kind=ikind)::nargin
    integer(kind=ikind)::i
    integer(kind=ikind),allocatable,dimension(:)::ipiv
    integer(kind=ikind), allocatable, dimension(:,:) :: limit_array
    integer(kind=ikind):: myid,numprocs,ierr
    integer(kind=ikind):: stat(mpi_status_size)
    character(len=200)::argin
    character(len=200)::fname,file_no
    character(len=200)::config_file ! character to hold file name
    character(len=200)::xsparse
    character(len=200)::method
    type(param)::prm
    !**********execution part**************************************************
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)
    nargin=iargc()
    if (nargin .ne. 1) then
       write(*,'(a)') 'usage: ./inverse_modeling.exe config_file'
       stop
    else
       call getarg (1,config_file)
    endif
    config_file=adjustl(trim(config_file))    
    prm=set_all_args(config_file)
    write(*,*)'*****************************************************'
    allocate(junk(1,1))
    junk(1,1)=0
    if (prm%covtype=='independent') then
       call read_binary_vector(prm%inp_direc,prm%filename_parameters,q_diag)
       allocate(q(prm%locations,sum(prm%days)*prm%withinday))
       q=reshape(q_diag,(/prm%locations,sum(prm%days)*prm%withinday/))
    else
       call read_binary_matrix(prm%inp_direc,prm%filename_parameters,parameters)
       prm%variance=parameters(:,1)
       prm%sclength=parameters(:,2)
       prm%tclengthd=parameters(:,3)
       allocate(temporal(sum(prm%days)*prm%withinday,sum(prm%days)*prm%withinday))
       temporal=constant_zero
       allocate(Q_sum(prm%size))
       call temporal_covariance(temporal,Q_sum,prm)
       call read_binary_matrix(prm%inp_direc,prm%filename_lat_long,sp_coord)
       allocate(hEff(prm%locations,prm%locations))
       if (myid .eq. 0) then
          hEff=distance_spherical(sp_coord,sp_coord,prm%locations,prm%locations)
       endif
       call mpi_bcast(hEff,prm%locations*prm%locations,mpi_real8,0,mpi_comm_world,ierr) 
    endif
    stime = mpi_wtime()
    allocate(temporary(prm%observations,prm%observations))       
    allocate (limit_array(numprocs,2))
    limit_array=limit(sum(prm%days)*prm%withinday,numprocs)
    temporary=0
    if (prm%covtype=='independent') then
       temporary=hqht_mkl_oc_diag(prm,myid,numprocs,limit_array(myid+1,1),limit_array(myid+1,2),q)
    else
       temporary=hqht_calc(prm%hincore,prm,myid,numprocs,limit_array(myid+1,1),limit_array(myid+1,2),temporal,hEff)
    endif
    call mpi_barrier(mpi_comm_world, ierr)
    if (myid .eq. 0) then
       allocate(hqht(prm%observations,prm%observations))
    endif
    call mpi_reduce(temporary,hqht,prm%observations*prm%observations,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    deallocate(limit_array,temporary)
    etime = mpi_wtime()
    if (myid .eq. 0) then
       fname='hqht.bin'
       call write_binary_matrix(prm%out_direc,fname,hqht) 
       write(*,*) 'total time taken for hqht: ', etime-stime , ': in seconds'
       fname=prm%filename_r
       call read_binary_vector(prm%inp_direc,fname,r)
       call psi(hqht,r)
	   !hqht=hqht+r        !hqht becomes hqhtr=psi
       fname='psi.bin'
       call write_binary_matrix(prm%out_direc,fname,hqht)
       deallocate(r)
    endif
    call mpi_finalize(ierr)
    end program inverse_modeling_main

