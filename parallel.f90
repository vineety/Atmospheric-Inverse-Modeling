    module parallel
    use constants
    use initialize
    use library_generic
    contains
    function limit(tcov_rows,numprocs) result (limit_array)
    ! purpose:
    ! load balancing function for hqht
    ! using sgemm blas routine
    !*******************variable definitions ********
    implicit none
    integer(kind=ikind),intent(in):: tcov_rows
    integer(kind=ikind),intent(in):: numprocs
    integer(kind=ikind):: remainder
    integer(kind=ikind):: divison
    integer(kind=ikind):: index1
    integer(kind=ikind):: pset
    integer(kind=ikind):: i
    integer(kind=ikind):: limit_array(numprocs,2)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    limit_array=0
    remainder=mod(tcov_rows,numprocs)
    if (remainder==0) then
        pset=0
    elseif (remainder>0) then
        pset=1
    endif
    divison=(tcov_rows-remainder)/numprocs
    index1=1
    do i = 1,numprocs
        if (remainder==0 .and. pset==0) then
            limit_array(i,1)= index1
            limit_array(i,2)= index1+divison-1
            index1=index1+divison
        elseif ((remainder>0 .or. remainder==0) .and. pset==1) then
            if (i>1) then
                index1=limit_array(i-1,2)+1
            endif
            if (remainder>0)then
                limit_array(i,1)=index1
                limit_array(i,2)=index1+divison-1+(remainder-(remainder-1))
                remainder=remainder-1
            elseif (remainder==0) then
                limit_array(i,1)=index1
                limit_array(i,2)=index1+divison-1
            endif
        endif
    enddo
    end function limit

    function month_query(strip_query,prm) result(month)
      integer(kind=ikind), intent(in) :: strip_query
      integer(kind=ikind) :: month
      type(param), intent(in) :: prm 
      integer(kind=ikind), allocatable,dimension (:,:) :: period
      allocate(period(prm%size,2))
      do i=1,prm%size
         if (i==1) then
            period(i,1)=1;
            period(i,2)=prm%days(i)*prm%withinday;
            if (strip_query>=period(i,1) .and. strip_query<= period(i,2)) then
               month = i;
               exit
            endif
         elseif (i>1) then
            period(i,1)=period(i-1,2)+1;
            period(i,2)=period(i-1,2)+(prm%days(i)*prm%withinday);
            if (strip_query>=period(i,1) .and. strip_query<= period(i,2)) then
               month = i;
               exit
            endif
         endif
      end do
    end function month_query

    end module parallel
