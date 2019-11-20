    module library_generic
    
    use constants
    
    contains
    
    function distance_spherical (x,y,xrows,yrows) result(distance)
    !**************************************************************************
    !written by : vineet yadav 
    !date : 8/29/2010
    !purpose : return spherical distance between two matrices having 2 cols
    !***********specification part*********************************************
    implicit none
    integer, parameter:: rkind8=kind(1.d0)
    integer(kind=ikind) :: i ! loop iterator
    integer(kind=ikind) :: xrows ! an integer(kind=ikind) specifying no of rows in y vector
    integer(kind=ikind):: yrows ! an integer(kind=ikind) specifying no of rows in y vector
    real(kind=rkind8) :: r = 6371.009 ! Mean radius of the earth International Union of Geodesy and Geophysica
    real(kind=rkind8) :: pi,lat1 ! value of pi
    real(kind=rkind8), dimension(:,:) :: distance(xrows,yrows) ! the matrix that would return the result
    real(kind=rkind8), dimension(:,:) :: x(xrows,2) ! the input matrix dimension; only two cols allowed
    real(kind=rkind8), dimension(:,:) :: y(yrows,2) ! the input matrix dimension; only two cols allowed
    real(kind=rkind8), dimension(:) :: d_temp(yrows) ! temporary variable
    real(kind=rkind8), dimension(:) :: dlat(yrows),dlat1(yrows) ! the input matrix dimension; only two cols allowed
    real(kind=rkind8), dimension(:) :: dlon(yrows),dlon1(yrows) ! the input matrix dimension; only two cols allowed
    real(kind=rkind8), dimension(:) :: lat2(yrows) ! the input matrix dimension; only two cols allowed
    !**********execution part**************************************************
    ! multiply x and y matrices by pi/180
    pi=3.14159265358979323846264338327950288419716939937510
    lat2=y(:,1)*pi/180
    !$omp parallel do private(i)schedule(static)
    do i=1,xrows
       dlat1 = x(i,1)-y(:,1)
       dlon1 = x(i,2)-y(:,2)
       dlat=dlat1*pi/180
       dlon=dlon1*pi/180
       lat1=x(i,1)*pi/180
       d_temp=(dsin(dlat/2))**2 + dcos(lat1)*dcos(lat2)*(dsin(dlon/2))**2
       distance(:,i)=r*2*atan2(sqrt(d_temp),sqrt(1-d_temp));
    enddo
    !$omp end parallel do
    end function distance_spherical
    
    function distance_euclidean (x,y,xrows, yrows) result(distance)
    !**************************************************************************
    !written by : vineet yadav 										   
    !date : 8/29/2010
    !purpose : returns euclidean distance between two matrices having 2 cols
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind) :: i ! loop iterator
    integer(kind=ikind) :: xrows ! an integer(kind=ikind) specifying no of rows in y vector
    integer(kind=ikind) :: yrows ! an integer(kind=ikind) specifying no of rows in y vector
    real(kind=rkind), dimension(:,:) :: distance(xrows,yrows) ! the matrix that would return the result
    real(kind=rkind), dimension(:,:) :: x(xrows,2) ! the input matrix dimension; only two cols allowed
    real(kind=rkind), dimension(:,:) :: y(yrows,2) ! the input matrix dimension; only two cols allowed
    !**********execution part**************************************************
    !$omp parallel do private(i)schedule(static)
    do i=1,xrows
        distance(i,:) = sqrt((x(i,1)-y(:,1))**2 + (x(i,2)-y(:,2))**2);
    enddo
    !$omp end parallel do
    end function distance_euclidean

    function sorted_array (sort_vector,rows) result(vector)
    !**************************************************************************
    !written by : vineet yadav  								   
    !date : 8/29/2010
    !purpose : sort a vector and return sorted vector by selection sort algo
    !***********specification part*********************************************
    integer(kind=ikind) :: rows ! number of rows in the vector
    integer(kind=ikind) :: i ! loop iterator
    integer(kind=ikind) :: j ! loop iterator
    integer(kind=ikind) :: l ! loop iterator
    integer(kind=ikind) :: c ! integer(kind=ikind) used by the algo
    real(kind=rkind) :: b
    real(kind=rkind), dimension(:):: vector(rows)
    real(kind=rkind), dimension(:):: sort_vector(rows)
    vector = sort_vector
    !$omp parallel do
    do i=l,rows
        b=vector(l)
        do j=l,rows
            if (b>=vector(j))then
                b=vector(j)
                c=j;
            endif
        enddo
        vector(c)=vector(l);
        vector(i)=b;
        l=l+1;
    enddo
    !$omp end parallel do
    end function sorted_array

    function tril (array,rows,cols,adiag) result (ltarray)
    !**************************************************************************
    !written by : vineet yadav								   
    !date : 8/29/2010
    !purpose : return lower triangular portion of a matrix
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind) :: rows
    integer(kind=ikind) :: cols
    integer(kind=ikind) :: i
    integer(kind=ikind) :: j
    integer(kind=ikind) :: adiag
    real(kind=rkind), dimension (:,:) :: array(rows,cols)
    real(kind=rkind), dimension (:,:) :: ltarray(rows,cols)
    !**********execution part**************************************************
    if (adiag>0 .or. adiag <= -rows) then
        print*, "adiag should be <=0 and > -rows"
    endif
    do i = 1,cols
        do j = i,rows
            if (j>=i-adiag) then
                ltarray(j,i) = array(j,i)
            endif
        enddo
    enddo
    end function tril

    function triu (array,rows,cols,adiag) result (utarray)
    !**************************************************************************
    !written by : vineet yadav								   
    !date : 8/29/2010
    !purpose : return lower triangular portion of a matrix
    !***********specification part*********************************************
    implicit none
    integer(kind=ikind) :: rows
    integer(kind=ikind) :: cols
    integer(kind=ikind) :: i
    integer(kind=ikind) :: j
    integer(kind=ikind) :: adiag
    real(kind=rkind), dimension (:,:) :: array(rows,cols)
    real(kind=rkind), dimension (:,:) :: utarray(rows,cols)
    !**********execution part**************************************************
    if (adiag<0 .or. adiag >= rows) then
        print*, "adiag should be >= 0 and < rows"
    endif
    do i = 1,cols
        do j = 1,i
            if (j<=i-adiag) then
                utarray(j,i) = array(j,i)
            endif
        enddo
    enddo
    end function triu

    function eye (rows) result(identity)
    !***************************************************************
    !written by : vineet yadav								   
    !date : 8/29/2010
    !purpose : create an identity matrix
    !***********specification part**********************************
    implicit none
    integer(kind=ikind) :: rows
    integer(kind=ikind) :: i
    integer(kind=ikind) :: j
    real(kind=rkind), dimension (:,:) :: identity(rows,rows)
    !**********execution part***************************************
    do i = 1,rows
        do j = 1,rows
            if (j==i) then
                identity(j,i) = 1
            else
                identity(j,i) = 0
            endif
        enddo
    enddo
    end function eye

    subroutine kron(k, a, b)
    !purpose :  
    !subroutine to compute kronecker product of two matrices
    !***********specification part**********************************
    implicit none
    integer(kind=ikind):: i
    integer(kind=ikind):: j
    integer(kind=ikind):: ma
    integer(kind=ikind):: na
    integer(kind=ikind):: mb
    integer(kind=ikind):: nb
    real(kind=rkind), intent(in):: a(:,:)
    real(kind=rkind), intent(in):: b(:,:)
    real(kind=rkind), intent(inout):: k(:,:)
    !**********execution part***************************************
    ma = ubound(a, 1)
    na = ubound(a, 2)
    mb = ubound(b, 1)
    nb = ubound(b, 2)
    if (size(k,1) /= ma*mb .or. size(k,2) /= na*nb) then
        write(*,*) 'k has invalid size'
    end if
    forall(i=1:ma, j=1:na)
        k(mb*(i-1)+1:mb*i,nb*(j-1)+1:nb*j) = a(i,j)*b
    end forall
    end subroutine kron
    
    function generate_numbers (range_a,range_b) result(c)
    !purpose :  
    !generate a vector from range_a to range_b
    !with increment=1
    !***********specification part**********************************
    implicit none
    integer(kind=ikind),intent(in)::range_a
    integer(kind=ikind),intent(in)::range_b
    integer(kind=ikind)::i
    real(kind=rkind),allocatable,dimension(:)::c
    allocate(c((range_b-range_a)+1))
    do i=range_a,range_b
        c(i)=i
    end do
    end function generate_numbers
    
    function zeros_vector(rows) result(vector)
    !purpose :  
    !create a all-zero vector with dimension as (rows)
    !***********specification part**********************************
    implicit none
    integer(kind=ikind)::rows
    real(kind=rkind),allocatable,dimension(:)::vector
    allocate(vector(rows))
    vector=0.0
    end function zeros_vector
    
    function zeros_matrix(rows,cols) result(matrix)
    !purpose :  
    !create a all-zero matrix with dimension as (rows,cols)
    !***********specification part**********************************
    implicit none
    integer(kind=ikind)::rows
    integer(kind=ikind)::cols
    real(kind=rkind),allocatable,dimension(:,:)::matrix
    allocate(matrix(rows,cols))
    matrix=0.0
    end function zeros_matrix
    
    function horz_concat(a,b) result(c)
    !purpose :  
    !horizontally concatenate matrix a & b
    !***********specification part**********************************
    implicit none
    real(kind=rkind),intent(in),dimension(:,:)::a
    real(kind=rkind),intent(in),dimension(:,:)::b
    real(kind=rkind),allocatable,dimension(:,:)::c
    allocate(c(size(a,1),size(a,2)+size(b,2)))
    c(1:size(a,1),1:size(a,2)) = a
    c(1:size(a,1),(size(a,2)+1):(size(a,2)+size(b,2))) = b
    end function horz_concat

    function vert_concat(a,b) result(c)
    !purpose :  
    !vertically concatenate matrix a & b
    !***********specification part**********************************
    implicit none
    real(kind=rkind),intent(in),dimension(:,:)::a
    real(kind=rkind),intent(in),dimension(:,:)::b
    real(kind=rkind),allocatable,dimension(:,:)::c
    allocate(c(size(a,1)+size(b,1),size(a,2)))
    c(1:size(a,1),1:size(a,2)) = a
    c((size(a,1)+1):(size(a,1)+size(b,1)),1:size(a,2)) = b
    end function vert_concat
    
    function vert_concat_vector(a,b) result(c)
    !purpose :  
    !vertically concatenate vector a & b
    !***********specification part**********************************
    implicit none
    real(kind=rkind),intent(in),dimension(:)::a
    real(kind=rkind),intent(in),dimension(:)::b
    real(kind=rkind),allocatable,dimension(:)::c
    allocate(c(size(a)+size(b)))
    c(1:size(a)) = a
    c((size(a)+1):(size(a)+size(b))) = b
    end function vert_concat_vector

    subroutine print_mat_2d(matrix,name_matrix)
    !***************************************************************
    ! written by : Mae Qiu
    ! date : 5/24/2013
    ! function to print 2d matrices
    ! variable assignment
    !***************************************************************
    character(len=10):: col
    character(len=20):: formt
    character(len=*), intent(in):: name_matrix
    integer(kind=ikind):: n
    integer(kind=ikind):: rows
    integer(kind=ikind):: columns
    real(kind=rkind), intent(in), dimension(:,:):: matrix
    ! execution part
    rows=size(matrix,1)    
    columns=size(matrix,2)

    write(col,"(i0)") columns
    col = trim(col)
    formt = '('//col//'f15.6))'
    write(*,"(a)"),' '
    write(*,"(a,a)"),'now printing: ', name_matrix
    write(*,"(a)"),' '
    do n = 1, rows
       write(*,formt) matrix(n,:)
    end do
    write(*,*),' '
    end subroutine print_mat_2d

    subroutine print_mat_1d(matrix,name_matrix)
    !***************************************************************
    ! written by : Mae Qiu
    ! date : 5/24/2013
    ! function to print 1d matrices
    ! variable assignment
    !***************************************************************
    character(len=20):: formt
    character(len=*), intent(in):: name_matrix
    integer(kind=ikind):: n
    integer(kind=ikind):: rows
    real(kind=rkind), intent(in), dimension(:):: matrix
    ! execution part

    rows=size(matrix,1)    
    formt = '(f15.6)'
    write(*,"(a)"),' '
    write(*,"(a,a)"),'now printing: ', name_matrix
    write(*,"(a)"),' '
    do n = 1, rows
       write(*,formt) matrix(n)
    end do
    write(*,*),' '
    end subroutine print_mat_1d
    
    end module library_generic
