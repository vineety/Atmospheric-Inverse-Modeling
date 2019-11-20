    module distance
    use constants
    
    !****************************************************************************
    !  modul: distance
    !
    !  purpose: three functions to compute spherical distance
    !****************************************************************************
    contains

    function dist_spher_cosines (lat,lon,rows) result(distance)
    !**************************************************************************
    !Written by : Mae Qiu & vineet yadav
    !Institution : Department of Global Ecology, Carnegie Institution for science
    !Date : 7/7/2013
    !Purpose : Return spherical/geodetic distance between two or more places in a 
    !matrix form by using the law of cosines based on spherical trigonometry. 
    !The polar radius and equitorial radius in this function are based on WGS 84 
    !model of the earth.
    !Reference: R. W. Sinnott, Virtues of the Haversine, Sky and Telescope, 
    !vol. 68, no. 2, 1984, p. 159
    !
    !Note: If you are using this function externally i.e., not as part of inverse 
    !modeling software then specify, ikind=4 in the main program and rkind = 4 or 8
    !4 in case of integer means 4 byte integers
    !4 and 8 in case of real numbers mean single or double precision real numbers
    !***********Specification of Variables*************************************
    implicit none
    integer(kind=ikind) :: i, j! loop iterator
    integer(kind=ikind) :: rows ! an integer(kind=ikind) specifying no of rows in y vector
    real(kind=rkind) :: equitorial_radius = 6378.1370 ! equitorial radius of the earth (WGS 84)
    real(kind=rkind) :: polar_radius = 6356.7523 ! polar radius of the earth (WGS 84)
    real(kind=rkind) :: pi = 3.14159265358979323846264338327950288 ! value of pi
    real(kind=rkind) :: mean_radius ! mean radius of the earth
    real(kind=rkind) :: radiance_const  
    real(kind=rkind), dimension(:) :: lat(rows) ! the input matrix dimension; only two cols allowed
    real(kind=rkind), dimension(:) :: lon(rows) ! the input matrix dimension; only two cols allowed
    real(kind=rkind), dimension(:,:) :: distance(rows,rows) ! the matrix that would return the result
    !**********execution part**************************************************
    ! for converting to radiance
    radiance_const=pi/180
    ! mean radius of the earth
    mean_radius = (2*equitorial_radius+polar_radius)/3 ! International Union of Geodesy and Geophysics formula
    ! multiply x and y matrices by pi/180
    lat = lat*radiance_const
    lon = lon*radiance_const
    ! compute sin and cosine of appropriate vectors
    do i=1,rows-1
        do j=i+1,rows
            distance(i,j) = mean_radius*acos(sin(lat(i))*sin(lat(j))+cos(lat(i)) &
            *cos(lat(j))*cos(lon(i)-lon(j)))
        enddo
    enddo
    distance=dble(distance)
    distance=distance+transpose(distance)
    end function dist_spher_cosines
    
    function dist_spher_haver (lat,lon,rows) result(distance)
    !**************************************************************************
    !Written by : vineet yadav
    !Institution : Department of Global Ecology, Carnegie Institution for science
    !Date : 7/7/2013
    !Purpose : Return spherical/geodetic distance between two or more places in a 
    !matrix form by using the haversine formula based on spherical trigonometry. 
    !The polar radius and equitorial radius in this function are based on WGS 84 
    !model of the earth.
    !Reference: R. W. Sinnott, Virtues of the Haversine, Sky and Telescope, 
    !vol. 68, no. 2, 1984, p. 159
    !
    !Note: If you are using this function externally i.e., not as part of inverse 
    !modeling software then specify, ikind=4 in the main program and rkind = 4 or 8
    !4 in case of integer means 4 byte integers
    !4 and 8 in case of real numbers mean single or double precision real numbers
    !***********Specification of Variables*************************************
    implicit none
    integer(kind=ikind) :: i, j ! loop iterators
    integer(kind=ikind) :: rows ! specifying no of rows in lat and long vectors
    real(kind=rkind) :: equitorial_radius = 6378.1370 ! equitorial radius of the earth (WGS 84)
    real(kind=rkind) :: polar_radius = 6356.7523 ! polar radius of the earth (WGS 84)
    real(kind=rkind) :: pi = 3.14159265358979323846264338327950288
    real(kind=rkind) :: mean_radius ! mean radius of the earth
    real(kind=rkind) :: radians_const ! radians conversion variable
    real(kind=rkind) :: lat1, dlat1, dlon1, dlat, dlon, d_temp ! temporary internal variables
    real(kind=rkind), dimension(:) :: lat(rows) ! the input matrix dimension; only one col is allowed
    real(kind=rkind), dimension(:) :: lon(rows) ! the input matrix dimension; only one col is allowed
    real(kind=rkind), dimension(:) :: lat2(rows)! temporary internal vector variable
    real(kind=rkind), dimension(:,:) :: distance(rows,rows) ! the matrix that would be returned as result
    !**********execution part**************************************************
    ! for converting to radians
    radians_const=pi/180
    ! mean radius of earth
    mean_radius = (2*equitorial_radius+polar_radius)/3  ! International Union of Geodesy and Geophysics formula
    ! multiply x and y matrices by pi/180 to convert to radians
    lat2=lat(:)*radians_const
    do i=1,rows-1
        do j= i+1,rows
            dlat1 = lat(i)-lat(j)
            dlon1 = lon(i)-lon(j)
            dlat=dlat1*radians_const
            dlon=dlon1*radians_const
            lat1=lat(i)*radians_const
            ! haversine formula
            d_temp=(sin(dlat/2))**2 + cos(lat1)*cos(lat2(j))*(sin(dlon/2))**2
            distance(i,j)=mean_radius*2*atan2(sqrt(d_temp),sqrt(1-d_temp));
        enddo
    enddo
    ! as only lower triangular matrix is computed we take the transpose of the lower triangular matrix
    ! and then add it to the lower triangular matrix to get the symmetrical distance matrix
    distance=distance+transpose(distance)
    end function dist_spher_haver

    function dist_spher_vincenty(lat,lon,rows) result (distance)
    !**************************************************************************
    !Written by : vineet yadav
    !Institution : Department of Global Ecology, Carnegie Institution for science
    !Date : 7/7/2013
    !Purpose : Return spherical/geodetic distance between two or more places in a 
    !matrix form by using the vincenty ellipsoidal iterative formula. 
    !The polar radius and equitorial radius in this function are based on WGS 84 
    !model of the earth.
    !Reference: T. Vincenty, Direct and Inverse Solutions of Geodesics on the Ellipsoid 
    !with application of nested equations, Survey Review. XXIII Vincenty, T. April 1975
    !
    !Note: If you are using this function externally i.e., not as part of inverse modeling
    !software then specify, ikind=4 in the main program and rkind = 4 or 8
    !4 in case of integer means 4 byte integers
    !4 and 8 in case of real numbers mean single or double precision real numbers
    !
    !Note: As distance is computed in meters in this function, in the end,
    !we divide by 1000 to convert it to Kms
    !
    !***********Specification of Variables*************************************
    implicit none
    integer(kind=ikind):: i, j , itercount ! loop iterators and counters
    integer(kind=ikind):: rows ! dimension of input vector lat and long
    integer(kind=ikind):: iterlimit=10! set limit for maximum iterations for vincenty formula
    real (kind=rkind) :: a = 6378137.0 ! in meters, equitorial axis
    real (kind=rkind) :: b = 6356752.31424518 ! in meters, polar axis wgs84 ellipsoid
    real (kind=rkind) :: f = 1/298.257223563 ! flattening coefficient
    real (kind=rkind) :: pi = 3.14159265358979323846264338327950288 ! value of pi
    real (kind=rkind) :: radians_const ! radians conversion variable
    real (kind=rkind) :: lambda, lambdaold, u1, u2, l, k2, lon1, lon2
    real (kind=rkind) :: c, sinsigma, cossigma, sigma, alpha ! temporary internal variables
    real (kind=rkind) :: temp_a, temp_b, deltasigma, cos2sigmam 
    real (kind=rkind), dimension(:) :: lat(rows) ! the input matrix dimension; only one col is allowed
    real (kind=rkind), dimension(:) :: lon(rows) ! the input matrix dimension; only one col is allowed
    real (kind=rkind), dimension(:,:) :: distance(rows,rows) ! the matrix that would be returned as result
    real (kind=rkind), dimension(:,:) :: dist(rows,rows) ! temporary internal matrix
    !**********execution part**************************************************
    ! convert inputs in degrees to radians:
    radians_const=pi/180
    distance=0.0
    lat = lat * radians_const
    lon = lon * radians_const
    do i=1,rows;
        ! correct for errors at exact poles by adjusting 0.6 millimeters:
        if (abs(pi/2-abs(lat(i))) < 1e-10) then
            lat(i) = signum_real(lat(i))*(pi/2-(1e-10))
        endif
    end do
    lat = atan((1-f)*tan(lat(:)));
    lon = mod(lon(:),2*pi);
    do i = 1,rows-1
        u1=lat(i)
        do j = i+1, rows
            u2=lat(j)
            l = abs(lon(i)-lon(j))
            if (l > pi) then
                l = 2*pi - l
            endif
            lambda=l
            lambdaold = 0
            itercount = 0
            do while (itercount /= iterlimit .or. abs(lambda-lambdaold) > 1e-12) 
                itercount = itercount+1
                if(itercount > iterlimit) then
                    !points are essentially antipodal. precision may be reduced slightly.
                    lambda = pi
                    exit
                endif
                lambdaold = lambda
                ! Vincenty's formula
                sinsigma = sqrt((cos(u2)*sin(lambda))**2+(cos(u1)* &
                sin(u2)-sin(u1)*cos(u2)*cos(lambda))**2)
                cossigma = sin(u1)*sin(u2)+cos(u1)*cos(u2)*cos(lambda)
                sinsigma=dble(sinsigma)
                cossigma=dble(cossigma)
                sigma = atan2(sinsigma,cossigma)
                alpha = asin(cos(u1)*cos(u2)*sin(lambda)/sin(sigma))
                cos2sigmam = cos(sigma)-2*sin(u1)*sin(u2)/cos(alpha)**2
                c = f/16*cos(alpha)**2*(4+f*(4-3*cos(alpha)**2))
                lambda = l+(1-c)*f*sin(alpha)*(sigma+c*sin(sigma)* &
                (cos2sigmam+c*cos(sigma)*(-1+2*cos2sigmam**2)))
                ! correct for convergence failure in the case of essentially antipodal points
                if(lambda > pi) then
                    !points are essentially antipodal. precision may be reduced slightly.
                    lambda = pi
                    exit
                endif
            end do
            k2 = cos(alpha)**2*(a**2-b**2)/b**2
            temp_a = 1+k2/16384*(4096+k2*(-768+k2*(320-175*k2)))
            temp_b = k2/1024*(256+k2*(-128+k2*(74-47*k2)))
            deltasigma = temp_b*sin(sigma)*(cos2sigmam+temp_b/4*(cos(sigma)*(-1+2*cos2sigmam**2) &
            -temp_b/6*cos2sigmam*(-3+4*sin(sigma)**2)*(-3+4*cos2sigmam**2)))
            dist(i,j) = b*temp_a*(sigma-deltasigma)
        end do
    end do
    ! as only lower triangular matrix is computed we take the transpose of the lower triangular matrix
    ! and then add it to the lower triangular matrix to get the symmetrical distance matrix
    distance=(dist+transpose(dist))/1000 ! convert distance to km
    end function dist_spher_vincenty

    function signum_real(a) result(b)
    !**************************************************************************
    !Written by : vineet yadav
    !Institution : Department of Global Ecology, Carnegie Institution for science
    !Date : 7/7/2013
    !Purpose : This is the standard signum mathematical function that returns -1 
    !if a real number is <0 and 0 if it is 0 and 1 if it is >1
    !
    !Note: If you are using this function externally i.e., not as part of inverse modeling
    !software then specify, ikind=4 in the main program and rkind = 4 or 8
    !4 in case of integer means 4 byte integers
    !4 and 8 in case of real numbers mean single or double precision real numbers
    !
    !***********Specification of Variables*************************************
    implicit none
    real(kind=rkind), intent(in) :: a ! input real number
    real(kind=rkind) :: b ! output value of signum function
    !**********execution part**************************************************
    if (a<0) then
        b=-1
    elseif (a==0) then
        b=0
    elseif (a>0) then
        b=1
    endif
    end function signum_real

    end module distance




