    module steps_mod
    !***********************************************************************
    !written by : Mae Qiu
    !date : 4/5/2013
    !purpose : wrapper functions to compute HQHt & HX, calling functions & 
    !          subroutines in library_inverse module
    !***********************************************************************
    use constants
    use library_generic
    use library_inverse
    use initialize

    contains
    
    function hqht_calc(in_core,prm,myid,numprocs,start,finish,temp_cov,hEff) result (hqht)
    ! purpose: compute hqht in core or out of core
    ! 
    !***********specification part*********************************************
    implicit none
    real(kind=rkind), intent(in), dimension(:,:)::temp_cov
    real(kind=rkind), intent(in), dimension(:,:)::hEff
    real(kind=rkind), allocatable, dimension(:,:)::sp_cov
    real(kind=rkind), allocatable, dimension(:,:)::sens_mat
    real(kind=rkind), allocatable, dimension(:,:)::hqht
    integer(kind=ikind),intent(in):: start
    integer(kind=ikind),intent(in):: finish
    integer(kind=ikind),intent(in):: myid,numprocs
    integer(kind=ikind):: tcov_rows
    integer(kind=ikind):: scov_rows
    integer(kind=ikind):: nnz
    character(len=*),intent(in)::in_core
    type(param),intent(in)::prm

    scov_rows=size(hEff,1)
    tcov_rows=size(temp_cov,1)

    allocate(hqht(prm%observations,prm%observations))
    allocate(sp_cov(size(hEff,1),size(hEff,2)))

    !***************compute hqht*******************************

    if (in_core=='yes') then
       !**********read in the full H matrix************************************
       allocate(sens_mat(prm%observations,tcov_rows*scov_rows))
       call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,0,prm%observations,tcov_rows*scov_rows,sens_mat,nnz)
       hqht=hqht_mkl_ic(prm,myid,numprocs,start,finish,temp_cov,hEff,sens_mat)
    elseif (in_core=='no') then
       hqht=hqht_mkl_oc_vc(prm,myid,numprocs,start,finish,temp_cov,hEff)
    else
       write(*,*) 'in_core is not defined properly, yes/no only!'
       stop
    endif

    end function hqht_calc

    function hx_calc(in_core,xsparse,prm,myid,numprocs,start,finish) result (hx)
    ! purpose: compute hx in core or out of core
    ! 
    !***********specification part*********************************************
    implicit none
    real(kind=rkind), allocatable, dimension(:,:)::h
    real(kind=rkind), allocatable, dimension(:,:)::x
    real(kind=rkind), allocatable, dimension(:,:)::temp
    real(kind=rkind), allocatable, dimension(:,:)::hx
    integer(kind=ikind),intent(in):: start
    integer(kind=ikind),intent(in):: finish
    integer(kind=ikind),intent(in):: myid,numprocs
    integer(kind=ikind):: nnz_h, nnz_x
    integer(kind=ikind):: x_rows
    character(len=*),intent(in)::in_core
    character(len=*),intent(in)::xsparse
    character(len=10)::fname
    type(param),intent(in)::prm

    allocate(hx(prm%observations,prm%covariates));
    x_rows=sum(prm%days)*prm%withinday*prm%locations
    !Both H & X are in core
    if (in_core=='yes') then       
        allocate(h(prm%observations,x_rows))
        allocate(x(x_rows,prm%covariates));
        call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,0,prm%observations,x_rows,h,nnz_h)
        !**********read in X based on it's sparse or not*************************
        if (xsparse=='yes') then
            call sparse_full_return(prm%inp_direc,prm%name_str_x,prm%ext,0,x_rows,prm%covariates,x,nnz_x)
        elseif (xsparse=='no') then
            fname=trim(adjustl(prm%name_str_x))//'f0'//trim(adjustl(prm%ext))
            call read_binary_matrix(prm%inp_direc,fname,temp)
            x=temp
            deallocate(temp)
        else
            write(*,*) 'xsparse is not defined properly, yes/no only!'
            stop
        endif        
        if ((nnz_h .ne. 0) .and. (nnz_x .ne. 0)) then
           hx=mat_mult(h,x)   !compute H*X
        endif
        deallocate(h,x)
    !Both H & X are out of core, read in H & X as strips
    elseif (in_core=='no') then
       hx=hx_atlas_oc(prm,myid,numprocs,start,finish,xsparse)
    else
       write(*,*) 'in_core is not defined properly, yes/no only!'
       stop
    endif

    end function hx_calc

    end module steps_mod
