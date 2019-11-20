    module library_inverse
    !****************************************************************************
    !  modul: library_inverse
    !  purpose: a series of functions and subroutines to do inverse modeling
    !****************************************************************************
    use constants
    use initialize
    use parallel
    use library_generic
    use file_io
    implicit none
    
    contains    
    subroutine psi(hqht,r)
    ! purpose:
    ! compute hqht+r, r is a 1-d array storing 
    ! diagonal elements for the original diagonal matrix r
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind)::i
    real(kind=rkind),intent(inout),dimension(:,:)::hqht
    real(kind=rkind),intent(in),dimension(:)::r
    !**********execution part**************************************************
    do i=1,size(r)
        hqht(i,i)=hqht(i,i)+r(i)
    end do
    end subroutine psi

    subroutine temporal_covariance(temporal_cov,Q_sum,prm)
    ! purpose:
    ! compute temporal covariance matrix
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind)::i
    real(kind=rkind),dimension(:)::Bound_1(2)
    real(kind=rkind),dimension(:)::Bound_2(2)
    real(kind=rkind),dimension(:)::Q_sum
    real(kind=rkind),dimension(:,:)::temporal_cov
    real(kind=rkind),allocatable,dimension(:,:)::temp_coord
    real(kind=rkind),allocatable,dimension(:,:)::temp_cov
    real(kind=rkind),allocatable,dimension(:,:)::tdist
    type(param),intent(in)::prm
    !**********execution part**************************************************
    if (prm%covtype .eq. 'spatial') then
        do i = 1,size(prm%days)
            if (i .eq. 1) then
                Q_sum(i)=(prm%days(i)*prm%withinday)-(prm%overlap*prm%withinday);
            elseif (i .gt. 1) then
                Q_sum(i)=prm%days(i)*prm%withinday
            end if
        end do
        temporal_cov = eye(sum(prm%days)*prm%withinday)
    elseif (prm%covtype .eq. 'space_time') then
        do i = 1,size(prm%days)
            allocate(temp_coord(prm%days(i),2))
            temp_coord(:,1)=generate_numbers(1,prm%days(i))
            temp_coord(:,2)=zeros_vector(prm%days(i))
            allocate(tdist(prm%days(i),prm%days(i)))
            tdist=distance_euclidean(temp_coord,temp_coord,prm%days(i),prm%days(i))
            tdist=exp(-tdist/prm%tclengthd(i));
            deallocate(temp_coord)
            allocate(temp_cov(prm%days(i)*prm%withinday,prm%days(i)*prm%withinday))
            call kron(temp_cov,tdist,eye(prm%withinday));
            deallocate(tdist)
            if (i .eq. 1) then
                Q_sum(i)=sum(temp_cov((prm%overlap*prm%withinday)+1:size(temp_cov,2),&
                (prm%overlap*prm%withinday)+1:size(temp_cov,2)));
                Bound_1(1)=1
                Bound_1(2)=size(temp_cov,1)
                Bound_2(1)=1
                Bound_2(2)=size(temp_cov,2)
            elseif (i .gt. 1) then
                Q_sum(i)=sum(temp_cov)
                Bound_1(1)=Bound_1(2)+1
                Bound_1(2)=Bound_1(2)+size(temp_cov,1)
                Bound_2(1)=Bound_2(2)+1
                Bound_2(2)=Bound_2(2)+size(temp_cov,2)
            endif
            temporal_cov(int(Bound_1(1)):int(Bound_1(2)),&
            int(Bound_2(1)):int(Bound_2(2)))=temp_cov
            deallocate(temp_cov)
        end do
    end if
    end subroutine temporal_covariance

    function hqht_mkl_ic(prm,myid,numprocs,start,finish,temp_cov,hEff,sens_mat) result (hqht)
    ! purpose:
    ! compute in core hqht
    ! calls sgemm/dgemm blas routine
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind),intent(in):: start
    integer(kind=ikind),intent(in):: finish
    integer(kind=ikind),intent(in):: myid,numprocs
    integer(kind=ikind):: counter,previous
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind):: tcov_rows
    integer(kind=ikind):: scov_rows
    integer(kind=ikind):: i,j
    real(kind=rkind), intent(in), dimension(:,:)::temp_cov
    real(kind=rkind), intent(in), dimension(:,:)::hEff
    real(kind=rkind), intent(in), dimension(:,:)::sens_mat
    real(kind=rkind), allocatable, dimension(:,:)::sp_cov
    real(kind=rkind), allocatable, dimension(:,:)::hqht
    real(kind=rkind), allocatable, dimension(:,:)::hq
    real(kind=rkind), allocatable, dimension(:,:)::senst
    real(kind=rkind), allocatable, dimension(:,:)::q
    type(param),intent(in)::prm
    !**********execution part**************************************************
    previous=0
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    scov_rows=size(hEff,1)
    tcov_rows=size(temp_cov,1)
    allocate(hqht(prm%observations,prm%observations))
    allocate(hq(prm%observations,scov_rows))
    allocate(senst(scov_rows,prm%observations))
    allocate(sp_cov(size(hEff,1),size(hEff,2)))
    hqht=0
    do i = start,finish
        counter=month_query(i,prm)
        if (counter /= previous) then
          sp_cov = prm%variance(counter)*exp(-hEff/prm%sclength(counter))
          previous=counter
        endif
        hq=0
        !$omp parallel do private(j) shared(sens_mat,temp_cov,i,hq)
        do j=1,tcov_rows
            if (temp_cov(j,i) .ne. 0 .and. temp_cov(j,i) .ne. 1) then
                hq=hq+(sens_mat(:,1+(scov_rows*j)-scov_rows:j*scov_rows)*temp_cov(j,i))
            elseif (temp_cov(j,i) .eq. 1) then
                hq=hq+sens_mat(:,1+(scov_rows*j)-scov_rows:j*scov_rows)
            endif
        enddo
        !$omp end parallel do
        hq=mat_mult(hq,sp_cov)
        senst=transpose(sens_mat(:,1+(scov_rows*i)-scov_rows:i*scov_rows))
        hqht=hqht+mat_mult(hq,senst)
        call system_clock(count=clock_end)
        write(*,900)'incore hqh_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        real(clock_end-clock_start)/real(clock_rate), ' seconds'
    enddo
    end function hqht_mkl_ic

    function hqht_mkl_oc_vc(prm,myid,numprocs,start,finish,temp_cov,hEff) result (hqht)
    ! purpose:
    ! compute out of core hqht
    ! calls sgemm/dgemm blas routine
    ! calls mkl_scoomm/mkl_dcoomm to perform sparse dense matrix multiplication
    !***********specification part*********************************************
    implicit none
    real(kind=rkind), allocatable, dimension (:,:):: hqht
    real(kind=rkind), allocatable, dimension (:,:):: sp_cov
    real(kind=rkind), allocatable, dimension (:,:):: hq
    real(kind=rkind), allocatable, dimension (:,:):: hqt
    real(kind=rkind), allocatable, dimension (:):: values
    !dec$ attributes align: 32 :: hqht,sp_cov,hq,hqt,values
    real(kind=rkind), intent(in), dimension (:,:):: temp_cov
    real(kind=rkind), intent(in), dimension (:,:):: hEff
    real(kind=rkind)::constant_one=1.0
    real(kind=rkind)::constant_zero=0.0
    integer(kind=ikind), allocatable, dimension(:,:):: coord
    integer(kind=ikind),intent(in):: start
    integer(kind=ikind),intent(in):: finish
    integer(kind=ikind),intent(in):: myid,numprocs
    integer(kind=ikind):: counter,previous
    integer(kind=ikind):: i
    integer(kind=ikind):: j
    integer(kind=ikind):: tcov_dim
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind)::nnz
    character(len=100)::fname
    character(len=100)::file_no
    type(param),intent(in)::prm
    !**********execution part**************************************************
    previous=0
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    tcov_dim=size(temp_cov,1)
    allocate(hqht(prm%observations,prm%observations))
    allocate(hq(prm%observations,prm%locations))
    allocate(hqt(prm%locations,prm%observations))
    allocate(sp_cov(size(hEff,1),size(hEff,2)))
    hqht=constant_zero
    !sp_cov = prm%variance(counter)*exp(-hEff/prm%sclength(counter))
    do i=start,finish
        counter=month_query(i,prm)
        if (counter /= previous) then
          sp_cov = prm%variance(counter)*exp(-hEff/prm%sclength(counter))
          previous=counter
        endif
        hq=constant_zero
        !$omp parallel do private(j) shared(temp_cov,i,hq)
        do j=1,tcov_dim
            if ((temp_cov(j,i) .ne. constant_zero) .and. (temp_cov(j,i) .ne. constant_one)) then
                call sparse_hq_return(prm%inp_direc,prm%name_str_h,prm%ext,j,hq,temp_cov(j,i))
            elseif (temp_cov(j,i) .eq. constant_one) then
                call sparse_hq_return(prm%inp_direc,prm%name_str_h,prm%ext,j,hq,constant_one)
            endif
        enddo
        !$omp end parallel do
        if (accuracy .eq. 'single') then
            call sgemm('n', 't', prm%locations, prm%observations, &
            prm%locations, constant_one, sp_cov, prm%locations, hq, &
            prm%observations, constant_zero, hqt, prm%locations)
        elseif(accuracy .eq. 'double') then
            call dgemm('n', 't', prm%locations, prm%observations, &
            prm%locations, constant_one, sp_cov, prm%locations, hq, &
            prm%observations, constant_zero, hqt, prm%locations)
        endif
        call coord_h(prm%inp_direc,prm%name_str_h,prm%ext,i,coord, values, nnz)

        if ((accuracy .eq. 'single') .and. (nnz .ne. 0)) then
            call mkl_scoomm('n', prm%observations, prm%observations, &
            prm%locations, constant_one, 'glnf', values, coord(1:nnz,1), &
            coord(1:nnz,2), nnz, hqt, prm%locations, constant_one, hqht, prm%observations)
            deallocate(values,coord)
        elseif((accuracy .eq. 'double') .and. (nnz .ne. 0)) then
            call mkl_dcoomm('n', prm%observations, prm%observations, &
            prm%locations, constant_one, 'glnf', values, coord(1:nnz,1), &
            coord(1:nnz,2), nnz, hqt, prm%locations, constant_one, hqht, prm%observations)
            deallocate(values,coord)
        endif

        call system_clock(count=clock_end)
        write(*,900)'hqh_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        real(clock_end-clock_start)/real(clock_rate), ' seconds'
    enddo
    deallocate(hqt)
    end function hqht_mkl_oc_vc

    function hqht_mkl_oc_diag(prm,myid,numprocs,start,finish,q) result (hqht)
    ! purpose:
    ! compute out of core hqht with diagonal Q
    ! calls sgemm/dgemm blas routine from within mat_mult2
    !***********specification part*********************************************
    implicit none
    real(kind=rkind), allocatable, dimension (:,:):: hqht
    real(kind=rkind), allocatable, dimension (:,:):: hq
    real(kind=rkind), allocatable, dimension (:,:):: h
    real(kind=rkind), intent(in), dimension (:,:):: q
    real(kind=rkind)::constant_zero=0.0
    integer(kind=ikind),intent(in):: start
    integer(kind=ikind),intent(in):: finish
    integer(kind=ikind),intent(in):: myid,numprocs
    integer(kind=ikind):: i
    integer(kind=ikind):: j
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind)::nnz
    character(len=100)::fname
    character(len=100)::file_no
    type(param),intent(in)::prm
    !**********execution part**************************************************
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    allocate(hqht(prm%observations,prm%observations))
    allocate(hq(prm%observations,prm%locations))
    allocate(h(prm%observations,prm%locations))
    hqht=constant_zero
    do i=start,finish
        hq=constant_zero
        call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,i,&
             prm%observations,prm%locations,h,prm%observations*prm%locations)
        do j=1,prm%locations
           hq(1:prm%observations,j)=h(1:prm%observations,j)*q(j,i)
        enddo
        if (i==1) then
           hqht=mat_mult2(hq,h,'n','t')
        else
           hqht=hqht+mat_mult2(hq,h,'n','t')
        endif

        call system_clock(count=clock_end)
        write(*,900)'hqh_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        real(clock_end-clock_start)/real(clock_rate), ' seconds'
    enddo
    write(*,*)'*****************************************************'
    deallocate(hq)
    end function hqht_mkl_oc_diag

    function hqstrip(temp_cov,hEff,prm,month,block_no) result (hqt)
    ! purpose:
    ! compute out of core hqt
    ! calls sgemm/dgemm blas routine
    ! calls mkl_scoomm/mkl_dcoomm to perform sparse dense matrix multiplication
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind),intent(in):: block_no
    integer(kind=ikind),intent(in)::month
    integer(kind=ikind):: j
    integer(kind=ikind):: tcov_dim
    real(kind=rkind)::constant_one=1.0
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind), intent(in), dimension (:,:):: temp_cov
    real(kind=rkind), intent(in),dimension(:,:)::hEff
    real(kind=rkind), allocatable,dimension(:,:)::sp_cov
    real(kind=rkind), allocatable, dimension (:,:):: hq
    real(kind=rkind), allocatable, dimension (:,:):: hqt
    type(param),intent(in)::prm
    !**********execution part**************************************************
    tcov_dim=size(temp_cov,1)
    allocate(hq(prm%observations,prm%locations),hqt(prm%observations,prm%locations))
    allocate(sp_cov(prm%locations,prm%locations))
    hq=0
    do j=1,tcov_dim
        if ((temp_cov(j,block_no) .ne. constant_zero) .and. (temp_cov(j,block_no) .ne. constant_one)) then
            call sparse_hq_return(prm%inp_direc,prm%name_str_h,prm%ext,j,hq,temp_cov(j,block_no))
        elseif (temp_cov(j,block_no) .eq. constant_one) then
            call sparse_hq_return(prm%inp_direc,prm%name_str_h,prm%ext,j,hq,constant_one)
        endif
    enddo
    sp_cov=prm%variance(month)*exp(-hEff/prm%sclength(month));
    if (accuracy .eq. 'single') then
        call sgemm('n', 'n', prm%observations, prm%locations, prm%locations, &
        constant_one, hq, prm%observations, sp_cov, prm%locations, &
        constant_zero, hqt, prm%observations)
    elseif(accuracy .eq. 'double') then
        call dgemm('n', 'n', prm%observations,prm%locations,prm%locations,&
        constant_one,hq,prm%observations,sp_cov,prm%locations,constant_zero&
        ,hqt,prm%observations)
    endif
    end function hqstrip

    function hqstrip_diag(q,prm,block_no) result (hq)
    ! purpose:
    ! compute out of core hq for diagonal Q
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind),intent(in):: block_no
    integer(kind=ikind):: j
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind), intent(in), dimension (:,:):: q
    real(kind=rkind), allocatable,dimension(:,:):: h
    real(kind=rkind), allocatable, dimension (:,:):: hq
    type(param),intent(in)::prm
    !**********execution part**************************************************
    allocate(h(prm%observations,prm%locations))
    allocate(hq(prm%observations,prm%locations))
    hq=constant_zero
    call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,block_no,&
             prm%observations,prm%locations,h,prm%observations*prm%locations)
    do j=1,prm%locations
       hq(:,j)=h(:,j)*q(j,block_no)
    enddo
    end function hqstrip_diag

    function h_x(prm,temp_cov,xsparse) result(hx)
    ! purpose:
    ! compute h*x out of core
    ! calls sgemm/dgemm blas routine
    !***********specification part*********************************************
    implicit none
    character(len=*),intent(in)::xsparse
    character(len=10)::fname
    character(len=20)::file_string
    integer(kind=ikind):: i
    integer(kind=ikind):: nnz_h
    integer(kind=ikind):: nnz_x
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind), allocatable, dimension(:,:):: coord
    real(kind=rkind)::constant_one=1.0
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind), allocatable, dimension (:):: values
    real(kind=rkind),intent(in),dimension(:,:)::temp_cov
    real(kind=rkind),allocatable,dimension(:,:)::hxt
    real(kind=rkind),allocatable,dimension(:,:)::h
    real(kind=rkind),allocatable,dimension(:,:)::hx
    real(kind=rkind),allocatable,dimension(:,:)::x
    real(kind=rkind),allocatable,dimension(:,:)::temp
    type(param),intent(in)::prm
    !******************imlement**********************
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    allocate(hxt(prm%covariates,prm%observations));
    hxt=0
    allocate(h(prm%observations,prm%locations));
    allocate(x(prm%covariates,prm%locations));
    allocate(hx(prm%observations,prm%covariates));
    hx=0
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    do i=1,size(temp_cov,1)
        x=0
        h=0
        call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,i,prm%observations,prm%locations,h,nnz_h)

        !**********read in X based on it's sparse or not*************************
        if (xsparse=='yes') then
            call sparse_full_return(prm%inp_direc,prm%name_str_x,prm%ext,i,prm%covariates,prm%locations,x,nnz_x)
        elseif (xsparse=='no') then
            write(file_string,'(i20)') i
            fname=trim(adjustl(prm%name_str_x))//'f'//trim(adjustl(file_string))//trim(adjustl(prm%ext))
            call read_binary_matrix(prm%inp_direc,fname,temp)
            x=temp
            deallocate(temp)
        else
            write(*,*) 'xsparse is not defined properly, yes/no only!'
            stop
        endif

        if ((nnz_h .ne. 0) .and. (nnz_x .ne. 0) .and. (accuracy .eq. 'single')) then
            call sgemm('n', 't', prm%covariates, prm%observations, &
            prm%locations, constant_one, x, prm%covariates, h, &
            prm%observations, constant_one, hxt, prm%covariates)
        elseif((nnz_h .ne. 0) .and. (nnz_x .ne. 0) .and. (accuracy .eq. 'double')) then
            call dgemm('n', 't', prm%covariates, prm%observations, &
            prm%locations, constant_one, x, prm%covariates, h, &
            prm%observations, constant_one, hxt, prm%covariates)
        endif
        call system_clock(count=clock_end)
        !write(*,900)'hx_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        !real(clock_end-clock_start)/real(clock_rate), ' seconds'
    end do
    write(*,*)'*****************************************************'
    hx=transpose(hxt)
    end function h_x

    function hx_atlas_oc(prm,myid,numprocs,start,finish,xsparse) result (hx)
    ! purpose:
    ! compute h*x with both h & x out of core
    ! calls sgemm/dgemm blas routine
    !***********specification part*********************************************
    implicit none
    character(len=*),intent(in)::xsparse
    character(len=10)::fname
    character(len=200)::file_string
    integer(kind=ikind),intent(in):: start
    integer(kind=ikind),intent(in):: finish
    integer(kind=ikind),intent(in):: myid,numprocs
    integer(kind=ikind):: i
    integer(kind=ikind):: nnz_h, nnz_x
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    real(kind=rkind), allocatable, dimension(:,:)::h
    real(kind=rkind), allocatable, dimension(:,:)::x
    real(kind=rkind), allocatable, dimension(:,:)::temp
    real(kind=rkind), allocatable, dimension(:,:)::hx
    type(param),intent(in)::prm
    !**********execution part**************************************************
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    allocate(h(prm%observations,prm%locations))
    allocate(x(prm%locations,prm%covariates))
    allocate(hx(prm%observations,prm%covariates))
    hx=0
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    do i = start,finish
        x=0
        h=0
        call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,i,prm%observations,prm%locations,h,nnz_h)
        !**********read in X based on it's sparse or not*************************
        if (xsparse=='yes') then
            call sparse_full_return(prm%inp_direc,prm%name_str_x,prm%ext,i,prm%locations,prm%covariates,x,nnz_x)
        elseif (xsparse=='no') then
            write(file_string,'(i20)') i
            fname=trim(adjustl(prm%name_str_x))//'f'//trim(adjustl(file_string))//trim(adjustl(prm%ext))
            call read_binary_matrix(prm%inp_direc,fname,temp)
            x=temp
            deallocate(temp)
        else
            write(*,*) 'xsparse is not defined properly, yes/no only!'
            stop
        endif

       !$omp parallel do private(j) shared(sens_mat,temp_cov,i,hq)
       if ((nnz_h .ne. 0) .and. (nnz_x .ne. 0)) then
          hx=hx+mat_mult(h,x)
       endif

       !$omp end parallel do
       call system_clock(count=clock_end)
       write(*,900)'hx_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
          real(clock_end-clock_start)/real(clock_rate), ' seconds'       
    enddo
    write(*,*)'*****************************************************'
    end function hx_atlas_oc

    function h_sp(prm,sp,temp_rows) result(hsp)
    ! purpose:
    ! compute h*sp for Bayesian method
    ! sp is the prior estimate of the flux distribution
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind):: i
    integer(kind=ikind):: nnz_h
    integer(kind=ikind):: nnz_s
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind),intent(in)::temp_rows
    integer(kind=ikind), allocatable, dimension(:,:):: coord
    real(kind=rkind)::constant_one=1.0
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind),allocatable, dimension (:):: values
    real(kind=rkind),allocatable,dimension(:,:)::h
    real(kind=rkind),allocatable,dimension(:)::hsp
    real(kind=rkind),allocatable,dimension(:)::s
    real(kind=rkind),intent(in),dimension(:)::sp
    type(param),intent(in)::prm
    !******************imlement**********************
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    allocate(h(prm%observations,prm%locations))
    allocate(s(prm%locations))
    allocate(hsp(prm%observations))
    hsp=0
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    do i=1,temp_rows
        h=0
        s=0
        s=sp(1+(i*prm%locations)-prm%locations:i*prm%locations)
        nnz_s=size(sp,1)
        call sparse_full_return(prm%inp_direc,prm%name_str_h,prm%ext,i,prm%observations,prm%locations,h,nnz_h)
        
        if ((nnz_h .ne. 0) .and. (nnz_s .ne. 0)) then
           hsp=hsp+vec_mult_single(h,s)
        endif

        call system_clock(count=clock_end)
        !write(*,900)'hsp_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        !real(clock_end-clock_start)/real(clock_rate), ' seconds'
    end do
    write(*,*)'*****************************************************'
    end function h_sp

    subroutine bmioc(prm,psi,hx,z,inv_bmi) 
    ! purpose:
    ! compute the inverse of big matrix
    ! inv_bmi =|(HQHt + R) H*X |-1
    !          |  (H*X)t    0  |  
    !***********specification part*********************************************
    implicit none
    character(len=500)::fname
    type(param),intent(in)::prm
    integer(kind=ikind)::info
    integer(kind=ikind)::nmax
    integer(kind=ikind)::n
    integer(kind=ikind)::p
    integer(kind=ikind):: i
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind),allocatable,dimension(:)::ipiv
    real(kind=rkind),intent(in),dimension(:)::z
    real(kind=rkind),intent(in),dimension(:,:)::psi
    real(kind=rkind),intent(in),dimension(:,:)::hx
    real(kind=rkind),allocatable,dimension(:,:)::temp
    real(kind=rkind),allocatable,dimension(:)::work
    real(kind=rkind),allocatable,dimension(:,:)::v_beta
    real(kind=rkind),allocatable,dimension(:)::beta
    real(kind=rkind),intent(inout),dimension(:,:)::inv_bmi
    real(kind=rkind),allocatable,dimension(:)::rhs_xsi
    
    !*******************store bmi for futher factorization**********************

    inv_bmi=vert_concat(horz_concat(psi,hx),horz_concat(transpose(hx),zeros_matrix(prm%covariates,prm%covariates))) 
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    if (accuracy .eq. 'single') then
        
        !************************calculations for inverse of v_beta******************
        call system_clock(count_rate=clock_rate) 
        call system_clock(count=clock_start)
        n=prm%observations
        p=prm%covariates
        allocate(ipiv(n),work(1))
        call ssytrf('u',n,psi,n,ipiv,work,-1,info)
        nmax=work(1)
        deallocate(work)
        allocate(work(nmax))
        call ssytrf('u',n,psi,n,ipiv,work,nmax,info)
        allocate(temp(n,p))
        temp=hx
        call ssytrs('u',n,p,psi,n,ipiv,temp,n,info)
        deallocate(ipiv,work)
        allocate(v_beta(p,p))
        v_beta = mat_mult(transpose(hx), temp)
        deallocate(temp)

        !************************Take Inverse of V_beta*******************************

        allocate(ipiv(p),work(1))
        call ssytrf('u',p,v_beta,p,ipiv,work,-1,info)
        nmax=work(1)
        deallocate(work)
        allocate(work(nmax))
        call ssytrf('u',p,v_beta,p,ipiv,work,nmax,info)
        allocate(temp(p,p))
        temp=eye(p)
        call ssytrs('u',p,p,v_beta,p,ipiv,temp,p,info)
        deallocate(ipiv,work)
        fname='v_beta.bin'
        call write_binary_matrix(prm%out_direc,fname,temp) 
        deallocate(v_beta,temp)

        !***************************Variance of beta computed*************************

        n=prm%observations+prm%covariates
        allocate(ipiv(n),work(1))
        allocate(rhs_xsi(n),beta(prm%covariates))
        call sgetrf(n,n,inv_bmi,n,ipiv,info)
        call sgetri(n,inv_bmi,n,ipiv,work,-1,info)
	    nmax=work(1)
	    deallocate(work)
	    allocate(work(nmax))
        call sgetri(n,inv_bmi,n,ipiv,work,nmax,info)
        rhs_xsi=vec_mult_single(inv_bmi,vert_concat_vector(z,zeros_vector(prm%covariates)))
        beta = rhs_xsi(prm%observations+1:n)
        fname='beta.bin'
        call write_binary_vector(prm%out_direc,fname,beta) 
        deallocate(beta,rhs_xsi,ipiv,work)
        call system_clock(count=clock_end)
        !write(*,900)'bmic ',1 ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        !real(clock_end-clock_start)/real(clock_rate), ' seconds'

    elseif (accuracy .eq. 'double') then
        
        call system_clock(count_rate=clock_rate) 
        call system_clock(count=clock_start)
        n=prm%observations
        p=prm%covariates
        allocate(ipiv(n),work(1))
        call dsytrf('u',n,psi,n,ipiv,work,-1,info)
        nmax=work(1)
        deallocate(work)
        allocate(work(nmax))
        call dsytrf('u',n,psi,n,ipiv,work,nmax,info)
        allocate(temp(n,p))
        temp=hx
        call dsytrs('u',n,p,psi,n,ipiv,temp,n,info)
        deallocate(ipiv,work)
        allocate(v_beta(p,p))
        v_beta = mat_mult(transpose(hx), temp)
        deallocate(temp)

        !************************Take Inverse of V_beta*******************************

        allocate(ipiv(p),work(1))
        call dsytrf('u',p,v_beta,p,ipiv,work,-1,info)
        nmax=work(1)
        deallocate(work)
        allocate(work(nmax))
        call dsytrf('u',p,v_beta,p,ipiv,work,nmax,info)
        allocate(temp(p,p))
        temp=eye(p)
        call dsytrs('u',p,p,v_beta,p,ipiv,temp,p,info)
        deallocate(ipiv,work)
        fname='v_beta.bin'
        call write_binary_matrix(prm%out_direc,fname,temp) 
        deallocate(v_beta,temp)

        !***************************Variance of beta computed*************************

        n=prm%observations+prm%covariates
        allocate(ipiv(n),work(1))
        allocate(rhs_xsi(n),beta(prm%covariates))
        call dgetrf(n,n,inv_bmi,n,ipiv,info)
        call dgetri(n,inv_bmi,n,ipiv,work,-1,info)
	    nmax=work(1)
	    deallocate(work)
	    allocate(work(nmax))
        call dgetri(n,inv_bmi,n,ipiv,work,nmax,info)
        rhs_xsi=vec_mult_single(inv_bmi,vert_concat_vector(z,zeros_vector(prm%covariates)))
        beta = rhs_xsi(prm%observations+1:n)
        fname='beta.bin'
        call write_binary_vector(prm%out_direc,fname,beta) 
        deallocate(beta,rhs_xsi,ipiv,work)
        call system_clock(count=clock_end)
        !write(*,900)'bmic ',1 ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        !real(clock_end-clock_start)/real(clock_rate), ' seconds'
    endif    
    write(*,*)'*****************************************************'    
    end subroutine bmioc

    subroutine qht_s(prm,temp_cov,heff,q,psi,inv_psi,hsp,z,qhts)
    ! purpose:
    ! compute qht*[inverse(hqht+r)*(z-h*sp)] out of core
    !***********specification part*********************************************
    implicit none
    character(len=500)::fname
    character(len=100)::file_no
    type(param),intent(in)::prm
    integer(kind=ikind)::info
    integer(kind=ikind)::temp_rows
    integer(kind=ikind)::period
    integer(kind=ikind)::counter
    integer(kind=ikind):: i
    integer(kind=ikind):: clock_start
    integer(kind=ikind):: clock_end
    integer(kind=ikind):: clock_rate
    integer(kind=ikind),allocatable,dimension(:)::ipiv
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind),intent(in),dimension(:)::z
    real(kind=rkind),intent(in),dimension(:)::hsp
    real(kind=rkind),intent(in),dimension(:,:)::psi
    real(kind=rkind),intent(in),dimension(:,:)::temp_cov
    real(kind=rkind),intent(in),dimension(:,:)::heff
    real(kind=rkind),intent(in),dimension(:,:)::q
    real(kind=rkind),intent(out),allocatable,dimension(:,:)::inv_psi
    real(kind=rkind),allocatable,dimension(:,:)::temp
    real(kind=rkind),allocatable,dimension(:)::zhsp
    real(kind=rkind),allocatable,dimension(:)::hqht_zhsp
    real(kind=rkind),allocatable,dimension(:)::work
    real(kind=rkind), allocatable, dimension (:,:):: hq
    real(kind=rkind), allocatable, dimension (:,:):: hqt
    real(kind=rkind), allocatable, dimension (:,:):: hqs
    real(kind=rkind),intent(out), allocatable, dimension (:):: qhts
    
    allocate(inv_psi(prm%observations,prm%observations))
    inv_psi=inverse_atlas(psi,'yes')
    allocate(zhsp(prm%observations))
    allocate(hqht_zhsp(prm%observations))
    zhsp=z-hsp       
    hqht_zhsp=vec_mult_single(inv_psi,zhsp)
    counter=1
    period = prm%days(counter)*prm%withinday
    temp_rows=sum(prm%days)*prm%withinday
    allocate(hq(prm%observations,prm%locations))
    allocate(hqt(prm%locations,prm%observations))
    allocate(hqs(prm%observations,prm%locations))
    allocate(qhts(prm%locations*temp_rows))
    hq=constant_zero
    hqt=constant_zero
    hqs=constant_zero
900 format (a,i6,a,i6,a,i6,a,f10.3,a)

    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    !****************DO QHt*[inverse(HQHt+R)*(z-H*Sprior)]********************
    do i=1,temp_rows ! compute for all strips
        if (prm%covtype=='independent') then
           hq=hqstrip_diag(q,prm,i)
        else
           hq=hqstrip(temp_cov,heff,prm,counter,i)
        endif
        hqt=transpose(hq)
        if (i .ge. (prm%overlap*prm%withinday)+1) then
           hqs=hqs+hq
           !hqts=hqts+hqt
        endif
        qhts(1+(i*prm%locations)-prm%locations:i*prm%locations)=&
                vec_mult_single(hqt,hqht_zhsp)
                !vec_mult_trans_single(hqt,hqht_zhsp,prm%observations)

        if (i .eq. period) then
           write(file_no,'(i20)') counter
           fname=adjustl(trim('hqsum_'//trim(adjustl(file_no))//'.bin'))
           call write_binary_matrix(prm%out_direc,fname,hqs) 
           counter=counter+1
           !if (counter .gt. size(prm%days)) then
           !   exit
           !endif
           hqs=constant_zero
           period=period+(prm%days(counter)*prm%withinday)
        endif

        call system_clock(count=clock_end)
        !write(*,900)'qhts_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        !        real(clock_end-clock_start)/real(clock_rate), ' seconds'
    enddo  
    write(*,*)'*****************************************************'
    end subroutine qht_s
   
    subroutine s_hat(prm,temp_cov,z,heff,q,inv_bmi,myid,numprocs,counter,&
                     start,finish,lambda_s,hqs,ms,xs,shat)
    ! purpose:
    ! compute shat for Geostatistical method
    !***********specification part*********************************************
    implicit none
    real(kind=rkind), allocatable, dimension (:,:):: lambda, lambda_s, hq, hqs
    real(kind=rkind), allocatable, dimension (:,:):: xo, x, xs, ms
    real(kind=rkind), dimension (:):: shat
    real(kind=rkind), intent(in), dimension(:,:)::heff, temp_cov, inv_bmi,q
    real(kind=rkind), intent(in), dimension(:):: z
    real(kind=rkind)::constant_one=1.0, constant_zero=0.0
    integer(kind=ikind), intent(in) :: myid, numprocs
    integer(kind=ikind), intent(in) :: counter, start, finish
    integer(kind=ikind):: i, period, tcov_dim
    integer(kind=ikind):: clock_start, clock_end, clock_rate
    integer(kind=ikind):: nnz,n,p
    integer(kind=ikind)::info
    type(param),intent(in)::prm
    character(len=100)::fname, file_no
    !**************************Execution******************
    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
900 format (a,i6,a,i6,a,i6,a,f10.3,a)
    tcov_dim=size(temp_cov,1)
    period = prm%days(counter)*prm%withinday
    ! for lambda summing
    allocate(lambda(prm%observations+prm%covariates,prm%locations))
    lambda=constant_zero
    lambda_s=constant_zero
     ! for hq summing
    allocate(hq(prm%observations,prm%locations))
    hqs=constant_zero
    hq=constant_zero
    ! for x summing
    allocate(xo(prm%locations,prm%covariates))
    allocate(x(prm%covariates,prm%locations))
    xo=constant_zero
    x=constant_zero
    xs=constant_zero
    ! for m summing
    ms=constant_zero
    ! for shat 
    shat=constant_zero
    n=prm%observations+prm%covariates
    p=prm%locations
    !*****************************All allocation completed******************************
    do i=start,finish ! compute for all strips
        if (prm%covtype=='independent') then
           hq=hqstrip_diag(q,prm,i)
        else
           hq=hqstrip(temp_cov,heff,prm,counter,i)
        endif
        call sparse_full_return(prm%inp_direc,prm%name_str_x,prm%ext,i,p,prm%covariates,xo,nnz)
        x=transpose(xo)  !transpose of the original x
        lambda=mat_mult(inv_bmi,vert_concat(hq,x))
        hqs=hqs+hq
        lambda_s=lambda_s+lambda(1:prm%observations,:)! this is for vshat
        ms=ms+lambda(prm%observations+1:size(lambda,1),:)
        xs=xs+x
        shat((1+(i*p)-p):(i*p))=vec_mult_single(transpose(lambda(1:prm%observations,:)),z)! create a s_hat vector
        call system_clock(count=clock_end)
        write(*,900)'shat_step ',i ,' done by node ',myid,' of ',numprocs,' time elapsed ',&
        real(clock_end-clock_start)/real(clock_rate), ' seconds'
    enddo
    end subroutine s_hat

    subroutine aposteriorig(prm,heff,q_sum)
    ! purpose:
    ! compute vshat and monthly/annual uncertainties
    ! for Geostatistical method
    !***********specification part*********************************************
    implicit none
    character(len=50)::fname
    character(len=50)::file_no
    integer(kind=ikind)::i
    integer(kind=ikind)::j
    integer(kind=ikind)::k
    integer(kind=ikind)::st
    integer(kind=ikind)::finish
    integer(kind=ikind)::clock_start
    integer(kind=ikind)::clock_end
    integer(kind=ikind)::clock_rate
    real(kind=rkind)::constant_one=1.0
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind),intent(in),dimension(:,:)::heff
    real(kind=rkind),intent(in),dimension(:)::q_sum
    real(kind=rkind),allocatable,dimension(:,:)::x
    real(kind=rkind),allocatable,dimension(:,:)::l_a
    real(kind=rkind),allocatable,dimension(:,:)::hq_a
    real(kind=rkind),allocatable,dimension(:,:)::m_a
    real(kind=rkind),allocatable,dimension(:,:)::x_a
    real(kind=rkind),allocatable,dimension(:,:)::q_a
    real(kind=rkind),allocatable,dimension(:,:)::xm
    real(kind=rkind),allocatable,dimension(:,:)::q_m
    real(kind=rkind),allocatable,dimension(:,:)::ms
    real(kind=rkind),allocatable,dimension(:,:)::hqs
    real(kind=rkind),allocatable,dimension(:,:)::lambda
    real(kind=rkind),allocatable,dimension(:,:)::month_uncert
    real(kind=rkind),allocatable,dimension(:,:)::ann_uncert
    type(param),intent(in)::prm
    !%%%%%%%%%%%%%%%%%%%execute
    allocate(l_a(prm%observations,prm%locations))
    allocate(hq_a(prm%observations,prm%locations))
    allocate(m_a(prm%covariates,prm%locations))
    allocate(x_a(prm%locations,prm%covariates))
    allocate(q_a(prm%locations,prm%locations))
    m_a=constant_zero
    l_a=constant_zero
    hq_a=constant_zero
    x_a=constant_zero
    q_a=constant_zero
900 format (a,i6,a,f10.3,a)

    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)
    do i = 1, size(prm%days)
        if (i .eq. 1) then
            st = (prm%overlap*prm%withinday)+1
            finish = prm%days(i)*prm%withinday
        elseif (i .gt. 1) then
            st =  finish+1
            finish = finish+prm%days(i)*prm%withinday
        endif

        write(file_no,'(i20)') i
        fname=adjustl(trim('m_'//trim(adjustl(file_no))//'.bin'))
        call read_binary_matrix(prm%out_direc,fname,ms)
        m_a = m_a+ms

        fname=adjustl(trim('x_'//trim(adjustl(file_no))//'.bin'))
        call read_binary_matrix(prm%out_direc,fname,xm)
        x_a = x_a+transpose(xm)

        fname=adjustl(trim('hqsum_'//trim(adjustl(file_no))//'.bin'))
        call read_binary_matrix(prm%out_direc,fname,hqs)
        hq_a=hq_a+hqs

        fname=adjustl(trim('lambda_'//trim(adjustl(file_no))//'.bin'))
        call read_binary_matrix(prm%out_direc,fname,lambda) 
        l_a=l_a+lambda

        allocate(q_m(prm%locations,prm%locations))
        q_m=constant_zero
        if (prm%covtype=='independent') then           
           do j=st,finish
              do k=1,prm%locations
                 q_m(k,k)=q_m(k,k)+q_sum((j-1)*prm%locations+k)
              enddo
           enddo
        else   
           q_m=(q_sum(i)*(prm%variance(i)*exp(-heff/prm%sclength(i))))
        endif
        q_a=q_a+q_m

        allocate(month_uncert(prm%locations,prm%locations))
        fname=adjustl(trim('month_uncert_'//trim(adjustl(file_no))//'.bin'))
        month_uncert=(mat_mult(-transpose(xm),ms)+q_m+&
        mat_mult(-transpose(hqs),lambda))/(finish-st+1)**2
        call write_binary_matrix(prm%out_direc,fname,month_uncert) 
        call system_clock(count=clock_end)
        write(*,900)'doing uncert ', i,' time elapsed ',&
                real(clock_end-clock_start)/real(clock_rate), ' seconds'
        deallocate(month_uncert,lambda,hqs,xm,ms,q_m)
    end do
    write(*,*)'*****************************************************'
    allocate(ann_uncert(prm%locations,prm%locations))
    ann_uncert=(mat_mult(-x_a,m_a)+q_a+mat_mult(-transpose(hq_a),l_a))&
    /((sum(prm%days)-prm%overlap)*prm%withinday)**2
    fname='ann_uncert.bin'
    call write_binary_matrix(prm%out_direc,fname,ann_uncert) 
    deallocate(ann_uncert)
    end subroutine aposteriorig

    subroutine aposteriorig_bayesian(prm,heff,q_sum,inv_psi)
    ! purpose:
    ! compute vshat and monthly/annual uncertainties
    ! for bayesian method
    !***********specification part*********************************************
    implicit none
    character(len=50)::fname
    character(len=50)::file_no
    integer(kind=ikind)::i
    integer(kind=ikind)::j
    integer(kind=ikind)::k
    integer(kind=ikind)::st
    integer(kind=ikind)::finish
    integer(kind=ikind)::clock_start
    integer(kind=ikind)::clock_end
    integer(kind=ikind)::clock_rate
    real(kind=rkind)::constant_one=1.0
    real(kind=rkind)::constant_zero=0.0
    real(kind=rkind),intent(in),dimension(:,:)::heff
    real(kind=rkind),intent(in),dimension(:)::q_sum
    real(kind=rkind),intent(in),dimension(:,:)::inv_psi
    real(kind=rkind),allocatable,dimension(:,:)::hq_a
    real(kind=rkind),allocatable,dimension(:,:)::q_a
    real(kind=rkind),allocatable,dimension(:,:)::q_m
    real(kind=rkind),allocatable,dimension(:,:)::hqs
    real(kind=rkind),allocatable,dimension(:,:)::month_uncert
    real(kind=rkind),allocatable,dimension(:,:)::ann_uncert
    type(param),intent(in)::prm
    !%%%%%%%%%%%%%%%%%%%execute
    allocate(hq_a(prm%observations,prm%locations))
    allocate(q_a(prm%locations,prm%locations))
    hq_a=constant_zero
    q_a=constant_zero
900 format (a,i6,a,f10.3,a)

    call system_clock(count_rate=clock_rate) 
    call system_clock(count=clock_start)    
    do i = 1, size(prm%days)
        if (i .eq. 1) then
            st = (prm%overlap*prm%withinday)+1
            finish = prm%days(i)*prm%withinday
        elseif (i .gt. 1) then
            st =  finish+1
            finish = finish+prm%days(i)*prm%withinday
        endif
        write(file_no,'(i20)') i
        fname=adjustl(trim('hqsum_'//trim(adjustl(file_no))//'.bin'))
        call read_binary_matrix(prm%out_direc,fname,hqs)
        hq_a=hq_a+hqs

        allocate(q_m(prm%locations,prm%locations))
        q_m=constant_zero
        if (prm%covtype=='independent') then           
           do j=st,finish
              do k=1,prm%locations
                 q_m(k,k)=q_m(k,k)+q_sum((j-1)*prm%locations+k)
              enddo
           enddo
        else   
           q_m=(q_sum(i)*(prm%variance(i)*exp(-heff/prm%sclength(i))))
        endif
        q_a=q_a+q_m

        allocate(month_uncert(prm%locations,prm%locations))
        fname=adjustl(trim('month_uncert_bayesian'//trim(adjustl(file_no))//'.bin'))
        month_uncert=(q_m-mat_mult(mat_mult(transpose(hqs),inv_psi),hqs))&
        /(finish-st+1)**2
        call write_binary_matrix(prm%out_direc,fname,month_uncert) 
        call system_clock(count=clock_end)
        write(*,900)'doing uncert ', i,' time elapsed ',&
                real(clock_end-clock_start)/real(clock_rate), ' seconds'
        deallocate(month_uncert,hqs,q_m)
    end do
    write(*,*)'*****************************************************'
    allocate(ann_uncert(prm%locations,prm%locations))
    ann_uncert=(q_a+mat_mult(mat_mult(-transpose(hq_a),inv_psi),hq_a))&
    /((sum(prm%days)-prm%overlap)*prm%withinday)**2
    fname='ann_uncert_bayesian.bin'
    call write_binary_matrix(prm%out_direc,fname,ann_uncert) 
    deallocate(ann_uncert)
    end subroutine aposteriorig_bayesian

    function vec_mult_trans_single(a,b,lda) result (c)
    ! a wrapper function to multiply two matrices using atlas or mkl
    ! atlas: automatically tuned linear algebra subroutines
    ! mkl: math kernel library
    ! as mkl is hardware specific library it would fastest in multiplying
    ! matrices
    ! to make this routine work you have to compile this fortran program
    ! by linking with atlas or mkl libraries		
    implicit none
    integer(kind=ikind):: a_rows ! no of rows in matrix a
    integer(kind=ikind):: a_cols ! no of cols in matrix a
    integer(kind=ikind):: b_cols=1 ! no of cols in matrix b
    integer(kind=ikind),intent(in):: lda ! leading dimension of matrix a
    real(kind=rkind),intent(in),dimension(:,:) :: a
    real(kind=rkind),intent(in),dimension(:) :: b
    real(kind=rkind),allocatable,dimension(:):: c
    real(kind=rkind):: alpha, beta ! two parameters required by sgemv subroutine 
    ! in atlas or mkl for multiplying matrices
    a_rows=size(a,1)
    a_cols=size(a,2)
    alpha=1.0
    beta=0.0
    allocate(c(a_rows))
    if (accuracy .eq. 'single') then
        call sgemv( 't', a_rows, a_cols, alpha, a, lda, b , 1, beta, c, 1  )
    elseif (accuracy .eq. 'double') then
        call dgemv( 't', a_rows, a_cols, alpha, a, lda, b , 1, beta, c, 1  )
    endif
    end function vec_mult_trans_single

    function vec_mult_single(a,b) result (c)
    ! a wrapper function to multiply two matrices using atlas or mkl
    ! atlas: automatically tuned linear algebra subroutines
    ! mkl: math kernel library
    ! as mkl is hardware specific library it would fastest in multiplying
    ! matrices
    ! to make this routine work you have to compile this fortran program
    ! by linking with atlas or mkl libraries		
    implicit none
    integer(kind=ikind):: a_rows ! no of rows in matrix a
    integer(kind=ikind):: a_cols ! no of cols in matrix a
    integer(kind=ikind):: b_cols=1 ! no of cols in matrix b
    real(kind=rkind),intent(in),dimension(:,:) :: a
    real(kind=rkind),intent(in),dimension(:) :: b
    real(kind=rkind),allocatable,dimension(:):: c
    real(kind=rkind):: alpha, beta ! two parameters required by sgemv subroutine 
    ! in atlas or mkl for multiplying matrices
    a_rows=size(a,1)
    a_cols=size(a,2)
    alpha=1.0
    beta=0.0
    allocate(c(a_rows))

    if (accuracy .eq. 'single') then
        call sgemv( 'n', a_rows, a_cols, alpha, a, a_rows, b , 1, beta, c, 1  )
    elseif (accuracy .eq. 'double') then
        call dgemv( 'n', a_rows, a_cols, alpha, a, a_rows, b , 1, beta, c, 1  )
    endif

    end function vec_mult_single

    function mat_mult(a,b) result (c)
    ! a wrapper function to multiply two matrices using atlas or mkl
    ! atlas: automatically tuned linear algebra subroutines
    ! mkl: math kernel library
    ! as mkl is hardware specific library it would fastest in multiplying
    ! matrices
    ! to make this routine work you have to compile this fortran program
    ! by linking with atlas or mkl libraries		
    implicit none
    integer(kind=ikind):: a_rows ! no of rows in matrix a
    integer(kind=ikind):: a_cols ! no of cols in matrix a
    integer(kind=ikind):: b_cols ! no of cols in matrix b
    real(kind=rkind),intent(in),dimension(:,:) :: a, b
    real(kind=rkind),allocatable,dimension(:,:):: c
    real(kind=rkind):: alpha, beta ! two parameters required by dgemm subroutine 
    ! in atlas or mkl for multiplying matrices
    a_rows=size(a,1)
    a_cols=size(a,2)
    b_cols=size(b,2)
    alpha=1.0
    beta=0.0
    allocate(c(a_rows,b_cols))
    if (accuracy .eq. 'single') then
        call sgemm('n', 'n', a_rows, b_cols, a_cols, alpha, a, &
        a_rows, b, a_cols, beta, c, a_rows)
    elseif(accuracy .eq. 'double') then
        call dgemm('n', 'n', a_rows, b_cols, a_cols, alpha, a, a_rows, b, a_cols, beta, c, a_rows)
    endif
    end function mat_mult

    function mat_mult2(a,b,transa,transb) result (c)
    ! a wrapper function to multiply two matrices using atlas or mkl
    ! atlas: automatically tuned linear algebra subroutines
    ! mkl: math kernel library
    ! as mkl is hardware specific library it would fastest in multiplying
    ! matrices
    ! to make this routine work you have to compile this fortran program
    ! by linking with atlas or mkl libraries            
    implicit none
    integer(kind=ikind):: a_rows ! no of rows in matrix a
    integer(kind=ikind):: a_cols ! no of cols in matrix a
    integer(kind=ikind):: b_rows ! no of rows in matrix b
    integer(kind=ikind):: b_cols ! no of cols in matrix b
    integer(kind=ikind):: m, n, k
    character(len=*),intent(in) :: transa, transb
    real(kind=rkind),intent(in),dimension(:,:) :: a, b
    real(kind=rkind),allocatable,dimension(:,:):: c
    real(kind=rkind):: alpha, beta ! two parameters required by dgemm subroutine 
    ! in atlas or mkl for multiplying matrices
    a_rows=size(a,1)
    a_cols=size(a,2)
    b_rows=size(b,1)
    b_cols=size(b,2)
    if (transa=='t') then
       m=a_cols
       k=a_rows
    else
       m=a_rows
       k=a_cols
    endif
    if (transb=='t') then
       n=b_rows
    else
       n=b_cols
    endif
    alpha=1.0
    beta=0.0
    allocate(c(m,n))
    if (rkind .eq. 4) then
        call sgemm(transa, transb, m, n, k, alpha, a, a_rows, b, b_rows, beta, c, m)
    elseif(rkind .eq. 8) then
        call dgemm(transa, transb, m, n, k, alpha, a, a_rows, b, b_rows, beta, c, m)
    endif
    end function mat_mult2

    function inverse_atlas(matrix,symmetric) result (matrix_inv)
    ! purpose:
    ! compute inverse of a square matrix
    ! calls ssytrf & ssytrs for symmetric matrix
    ! calls sgetrf & sgetri for non-symmetric matrix
    !***********specification part*********************************************
    implicit none
    real(kind=rkind), intent(in), dimension(:,:)::matrix
    real(kind=rkind), allocatable, dimension(:,:)::matrix_inv
    real(kind=rkind), allocatable, dimension(:)::work
    real(kind=rkind), allocatable, dimension(:,:)::temp
    real(kind=rkind), dimension(:)::qwork(1)  !workspace query to find out recommended lwork value
    integer(kind=ikind),allocatable, dimension(:):: ipiv
    integer(kind=ikind),allocatable, dimension(:):: gipiv
    integer(kind=ikind):: info
    integer(kind=ikind):: lwork
    integer(kind=ikind):: rows
    character(len=*),intent(in)::symmetric
    
    rows=size(matrix,1)
    allocate (matrix_inv(rows,rows))
    allocate (ipiv(rows)) 
    allocate (gipiv(rows)) 
    allocate (temp(rows,rows)) 

    !**********Compute the inverse of the matrix*****************************

    !Step1: call ?getrf to compute the LU factorization of a general m-by-n matrix A as A = P*L*U
    !       or call ?sytrf for sysmetric matrix 
    !Step2: call ?getri to compute the inverse of an LU-factored general matrix A (A has been
    !       overwritten after step1)
    if (rkind .eq. 4) then
       if (symmetric=='yes') then
        !***step1:
          call ssytrf('u',rows,matrix,rows,gipiv,qwork,-1,info)
          lwork=qwork(1)
          allocate(work(lwork))
          call ssytrf('u',rows,matrix,rows,ipiv,work,lwork,info)
          deallocate(work)
        !***step2:
          temp=eye(rows)
          call ssytrs('u',rows,rows,matrix,rows,ipiv,temp,rows,info)
          matrix_inv=temp
       elseif (symmetric=='no') then
        !***step1:
          call sgetrf(rows,rows,matrix,rows,ipiv,info)
        !***step2:
          call sgetri(rows,matrix,rows,gipiv,qwork,-1,info)  ! workspace query 
          lwork=qwork(1)
          allocate(work(lwork))
          call sgetri(rows,matrix,rows,ipiv,work,lwork,info)
          deallocate(work)          
          matrix_inv=matrix
       else
          write(*,*) 'symmetric is not defined properly, yes/no only!'
          stop
       endif 
    elseif(rkind .eq. 8) then
       if (symmetric=='yes') then
        !***step1:
          call dsytrf('u',rows,matrix,rows,gipiv,qwork,-1,info)
          lwork=qwork(1)
          allocate(work(lwork))
          call dsytrf('u',rows,matrix,rows,ipiv,work,lwork,info)
          deallocate(work)
        !***step2:
          temp=eye(rows)
          call dsytrs('u',rows,rows,matrix,rows,ipiv,temp,rows,info)
          matrix_inv=temp
       elseif (symmetric=='no') then
        !***step1:
          call dgetrf(rows,rows,matrix,rows,ipiv,info)
        !***step2:
          call dgetri(rows,matrix,rows,gipiv,qwork,-1,info)  ! workspace query 
          lwork=qwork(1)
          allocate(work(lwork))
          call dgetri(rows,matrix,rows,ipiv,work,lwork,info)
          deallocate(work)          
          matrix_inv=matrix
       else
          write(*,*) 'symmetric is not defined properly, yes/no only!'
          stop
       endif
    endif
    deallocate(gipiv,ipiv,temp)
    end function inverse_atlas

    end module library_inverse
