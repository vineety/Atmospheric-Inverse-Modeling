    module file_io
    use constants

    implicit none

    contains

    subroutine read_binary_matrix(inp_direc,file_name,matrix) 
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : read 2-d array written in binary format
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: fname 
    integer(kind=ikind):: rows     !number of rows 
    integer(kind=ikind):: cols     !number of cols
    integer(kind=ikind):: iocheck
    integer(kind=ikind),parameter :: unitcons=200
    real(kind=rkind),allocatable, dimension(:,:)::matrix
    !**********execution part**************************************************
    fname=trim(adjustl(inp_direc))//trim(adjustl(file_name))
    open (unitcons,file=fname,form='unformatted',access='stream',&
          iostat=iocheck,status="old",action="read",err=100)
    read(unitcons,iostat=iocheck,err=100) rows
    read(unitcons,iostat=iocheck,err=100) cols
    allocate(matrix(rows,cols))
    read(unitcons,iostat=iocheck,err=100) matrix
    close(unitcons)
100 call error_check(iocheck,fname)
    end subroutine read_binary_matrix
    
    subroutine read_binary_vector(inp_direc,file_name,vector) 
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : read 1-d vector written in binary format
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: fname 
    integer(kind=ikind):: elements ! number of elements for the vector
    integer(kind=ikind):: iocheck
    integer(kind=ikind),parameter :: unitcons=200
    real(kind=rkind),allocatable, dimension(:)::vector
    !**********execution part**************************************************
    fname=trim(adjustl(inp_direc))//trim(adjustl(file_name))!file name with directory
    open (unitcons,file=fname,form='unformatted',access='stream',&
          iostat=iocheck,status="old",action="read",err=100)
    read(unitcons,iostat=iocheck,err=100) elements
    allocate(vector(elements))
    read(unitcons,iostat=iocheck,err=100) vector
    close(unitcons)
100 call error_check(iocheck,fname)
    end subroutine read_binary_vector
    
    subroutine read_rows_cols(inp_direc,file_name,rows,cols) 
    !**************************************************************************
    !written by : mae qiu
    !date : 8/5/2013
    !purpose : read the rows & cols from the binary data file
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: fname 
    integer(kind=ikind):: rows     !number of rows 
    integer(kind=ikind):: cols     !number of cols
    integer(kind=ikind):: iocheck
    integer(kind=ikind),parameter :: unitcons=200
    !**********execution part**************************************************
    fname=trim(adjustl(inp_direc))//trim(adjustl(file_name))
    open (unitcons,file=fname,form='unformatted',access='stream',&
          iostat=iocheck,status="old",action="read",err=100)
    read(unitcons,iostat=iocheck,err=100) rows
    read(unitcons,iostat=iocheck,err=100) cols
    close(unitcons)
100 call error_check(iocheck,fname)
    end subroutine read_rows_cols
    
    subroutine write_binary_matrix(inp_direc,file_name,matrix) 
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : write 2-d array in binary format
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: name
    integer (kind=ikind):: rows 
    integer (kind=ikind):: cols
    integer (kind=ikind),parameter :: unitcons=200
    real (kind=rkind),intent(in), dimension(:,:)::matrix
    !**********execution part**************************************************
    rows=size(matrix,1)
    cols=size(matrix,2)
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))!file name with directory
    open (unitcons,file=name,form='unformatted',access='stream')
    write(unitcons) rows, cols, matrix
    close(unitcons)
    end subroutine write_binary_matrix

    subroutine write_binary_vector(inp_direc,file_name,vector) 
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : write 1-d vector in binary format
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: name
    integer (kind=ikind):: rows 
    integer (kind=ikind),parameter :: unitcons=200
    real (kind=rkind),intent(in), dimension(:)::vector
    !**********execution part**************************************************
    rows=size(vector,1)
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))!file name with directory
    open (unitcons,file=name,form='unformatted',access='stream')
    write(unitcons) rows, vector
    close(unitcons)
    end subroutine write_binary_vector
    
    subroutine write_binary_vector_int(inp_direc,file_name,vector) 
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : write 1-d integer vector in binary format
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: name
    integer (kind=ikind):: rows !number of cols in the file for sparse matrix
    integer (kind=ikind),parameter :: unitcons=200
    integer (kind=rkind),intent(in), dimension(:)::vector
    !**********execution part**************************************************
    rows=size(vector,1)
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))!file name with directory
    open (unitcons,file=name,form='unformatted',access='stream')
    write(unitcons) rows, vector
    close(unitcons)
    end subroutine write_binary_vector_int

    subroutine file_return_array(inp_direc,file_name,rows,cols,skip_lines,array) 
    !**************************************************************************
    !written by :vineet yadav  									   
    !date : 8/29/2010
    !purpose : read csv or tab delimited file and return an array
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: name
    integer(kind=ikind) :: rows !number of lines in the file
    integer(kind=ikind) :: cols !number of cols in the file
    integer(kind=ikind) :: skip_lines! number of lines you want to skip in the file
    integer(kind=ikind) :: i ! loop iterator
    real(kind=rkind),intent(inout),dimension(:,:) :: array(rows-skip_lines,cols)
    !**********execution part**************************************************
    if (rows < 1 .or. cols < 1 .or. skip_lines <0 ) then
        print*,"rows, cols and skip_lines have to be positive"
        stop
    endif
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))
    ! open the specified file
    open(unit=1,file=name,form="formatted",status="old",action="read")
    ! skip_lines
    if (skip_lines>0) then
        do i=1,skip_lines
            read(1,*) 
        enddo
    elseif (skip_lines == 0) then
        ! read file in an array
        do i=1,rows-skip_lines
            read(1,*) array(i,:)
        enddo
    endif
    close(1)
    end subroutine file_return_array

    subroutine sparse_write_coo(out_direc,file_name,nnz,coord,values)
    !**************************************************************************
    !written by :vineet yadav
    !date : 8/29/2010
    !purpose : write a binary coordinate file from sparse coordinates
    !***********specification part*********************************************
    implicit none
    character (len=*),intent(in):: out_direc ! output directory for writing file
    character (len=*),intent(in):: file_name ! output file name
    character (len=200):: name ! combined name of out_direc and file_name
    integer(kind=2),parameter :: unitcons=200
    integer(kind=ikind),intent(in) :: nnz ! total non zero sparse entries
    integer(kind=ikind),intent(in),dimension(:,:)::coord(nnz,2)
    real(kind=rkind),intent(in),dimension(:):: values(nnz)
    !**********execution part**************************************************
    if (size(coord,2) < 2 .or. size(coord,2) >2 ) then
        print*,"only two columns are allowed for the locations of non zero values in a matrix"
        stop
    endif
    name=trim(adjustl(out_direc))//trim(adjustl(file_name))
    open (unitcons,file=name,form='unformatted',access='stream')
    write(unitcons) nnz
    write(unitcons) coord,values
    close(unitcons)
    end subroutine sparse_write_coo

    subroutine sparse_read_coo(inp_direc,file_name,nnz,coord,values)
    !**************************************************************************
    !written by :vineet yadav
    !date : 8/29/2010
    !purpose : read a binary coordinate file from sparse coordinates
    !***********specification part*********************************************
    implicit none
    character (len=*),intent(in):: inp_direc ! output directory for writing file
    character (len=*),intent(in):: file_name ! output file name
    character (len=200):: name ! combined name of out_direc and file_name
    integer(kind=2),parameter :: unitcons=200
    integer(kind=2),parameter :: cols=2
    integer(kind=ikind) :: nnz ! total non zero sparse entries
    integer(kind=ikind),allocatable, dimension(:,:)::coord
    real(kind=rkind),allocatable,dimension(:):: values
    !**********execution part**************************************************
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))
    open (unitcons,file=name,form='unformatted',access='stream')
    read(unitcons) nnz
    allocate(values(nnz),coord(nnz,cols))
    read(unitcons) coord, values
    close(unitcons)
    end subroutine sparse_read_coo
    
    subroutine coord_h(inp_direc,name_str,ext,unitno,coord,values,nnz)
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : read in data from a coordinate format sparse file
    !***********specification part*********************************************
    implicit none
    character (len=*),intent(in):: inp_direc
    character (len=*),intent(in):: name_str
    character (len=*),intent(in):: ext
    character(len=200)::fname
    character(len=200)::file_string
    integer(kind=ikind),parameter :: cols=2 !number of cols in the file for sparse matrix
    integer(kind=ikind),parameter :: unitcons=200
    integer(kind=ikind),intent(in):: unitno ! unit no for opening files
    integer(kind=ikind)::i
    integer(kind=ikind)::row_index
    integer(kind=ikind):: iocheck
    integer(kind=ikind),intent(inout)::nnz
    integer(kind=ikind),intent(inout),allocatable,dimension(:,:)::coord
    real(kind=rkind),intent(inout),allocatable,dimension(:)::values
    !**********execution part**************************************************
    write(file_string,'(i20)') unitno
    fname=trim(adjustl(inp_direc))//trim(adjustl(name_str))//&
    trim(adjustl(file_string))//trim(adjustl(ext))
    open (unitcons+unitno,file=fname,form='unformatted',access='stream',&
          iostat=iocheck,status="old",action="read",err=100)
    read(unitcons+unitno,iostat=iocheck,err=100) nnz
    if (nnz .ne. 0) then
        allocate(values(nnz),coord(nnz,cols))
        read(unitcons+unitno,iostat=iocheck,err=100) coord, values
        close(unitcons+unitno)
        if (values(nnz) .eq. 0) then
            nnz=nnz-1
        endif
    else
       close(unitcons+unitno)
    endif
100 call error_check(iocheck,fname)
    end subroutine coord_h
    
    subroutine sparse_hq_return(inp_direc,name_str,ext,unitno,hq,multiplier) 
    !**************************************************************************
    !written by :vineet yadav
    !date : 8/29/2010
    !purpose : return coordinate from sparse array from a binary file
    !***********specification part*********************************************
    implicit none
    character (len=*),intent(in):: inp_direc
    character (len=*),intent(in):: name_str
    character (len=*),intent(in):: ext
    character(len=200)::fname
    character(len=200)::file_string
    integer(kind=ikind):: nnz
    integer(kind=ikind):: i !number of rows in the file or
    integer(kind=ikind):: iocheck
    integer(kind=ikind),allocatable,dimension(:,:)::coord
    integer(kind=2),parameter :: cols=3 !number of cols in the file
    integer(kind=2),parameter :: unitcons=200
    integer(kind=ikind),intent(in)::unitno
    real(kind=rkind),intent(in)::multiplier
    real(kind=rkind),intent(inout),dimension(:,:) :: hq
    real(kind=rkind),allocatable,dimension(:) :: values
    !**********execution part**************************************************
    write(file_string,'(i20)') unitno
    fname=trim(adjustl(inp_direc))//trim(adjustl(name_str))//&
    trim(adjustl(file_string))//trim(adjustl(ext))
    open (unitcons+unitno,file=fname,form='unformatted',access='stream',&
          iostat=iocheck,status="old",action="read",err=100)
    read(unitcons+unitno,iostat=iocheck,err=100) nnz
    if (nnz .ne. 0) then
        allocate(values(nnz),coord(nnz,cols-1))
        read(unitcons+unitno,iostat=iocheck,err=100) coord, values
        close(unitcons+unitno)
        if (multiplier .ne. 0 .and. multiplier  .ne. 1) then
            values=values*multiplier
            do i = 1, nnz
                hq(coord(i,1),coord(i,2))= hq(coord(i,1),coord(i,2))+values(i)
            enddo
        elseif (multiplier .eq. 1) then
            do i = 1, nnz
                hq(coord(i,1),coord(i,2))=hq(coord(i,1),coord(i,2))+values(i)
            enddo
        endif
    else
        close(unitcons+unitno)    
    endif
100 call error_check(iocheck,fname)
    end subroutine sparse_hq_return

    subroutine sparse_full_return(inp_direc,name_str,ext,unitno,dim1,dim2,matrix,nnz) 
    !**************************************************************************
    !written by : vineet yadav
    !date : 8/29/2010
    !purpose : use coordinate format sparse file to construct full matrix
    !***********specification part*********************************************
    implicit none
    character (len=*),intent(in):: inp_direc
    character (len=*),intent(in):: name_str
    character (len=*),intent(in):: ext
    character(len=200)::fname
    character(len=200)::file_string
    integer(kind=2),parameter :: unitcons=200
    integer(kind=ikind),intent(in):: unitno ! unit no for opening files
    integer(kind=ikind),intent(in):: dim1
    integer(kind=ikind),intent(in):: dim2
    integer(kind=ikind)::i
    integer(kind=ikind)::row_index
    integer(kind=ikind)::nnz
    integer(kind=ikind):: iocheck
    integer(kind=ikind),allocatable::coord(:,:)
    real(kind=rkind),allocatable::values(:)
    real(kind=rkind)::matrix(dim1,dim2)
    !**********execution part**************************************************
    write(file_string,'(i20)') unitno
    fname=trim(adjustl(inp_direc))//trim(adjustl(name_str))//&
    trim(adjustl(file_string))//trim(adjustl(ext))
    open (unitcons+unitno,file=fname,form='unformatted',access='stream',&
          iostat=iocheck,status="old",action="read",err=100)
    read(unitcons+unitno,iostat=iocheck,err=100) nnz
    if (nnz .ne. 0) then
        allocate(values(nnz),coord(nnz,2))
        read(unitcons+unitno,iostat=iocheck,err=100) coord, values
        close(unitcons+unitno)
        matrix=0
        row_index=size(coord,1)
        !$omp parallel do schedule(static) private(i) shared(matrix,values)
        do i = 1, row_index
            matrix(coord(i,1),coord(i,2))=values(i);
        enddo
        !$omp end parallel do
    elseif (nnz .eq. 0) then
        matrix=0
        close(unitcons+unitno)   
    endif
100 call error_check(iocheck,fname)
    end subroutine sparse_full_return

    function total_lines(inp_direc, file_name) result (total_lines_in_file)
    !**************************************************************************
    !written by : vineet yadav  									   
    !date : 8/29/2010
    !purpose : count number of rows in a comma or tab seperated file
    !***********specification part*********************************************
    implicit none
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: name
    integer(kind=ikind) :: totl_lines = 0 !initialize total number of lines
    integer(kind=ikind) :: total_lines_in_file ! variable that would be returned 
    integer(kind=ikind) :: i ! loop iterator
    !**********execution part**************************************************
    ! open the specified file
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))!file name with directory
    open(unit=1,file=name,form="formatted",status="old",action="read")
    ! count number of rows in the formatted file
    do
        read(1,*,iostat=i) 
        if (i==0) then
            totl_lines=totl_lines+1
        elseif (i/=0) then
            exit
        endif
    enddo
    close(1)
    total_lines_in_file = totl_lines
    end function total_lines

    function total_columns(inp_direc, file_name) result(total_columns_in_file)
    !**************************************************************************
    !written by : vineet yadav  										   
    !date : 8/29/2010
    !purpose : count number of columns in a comma, tab or space seperated file
    !***********specification part*********************************************
    implicit none
    character :: buffer !  character variable to hold a character
    character (len=*) :: inp_direc ! input directory from which to read file
    character (len=*) :: file_name ! name of the file to be read
    character (len=200) :: name
    integer(kind=ikind) :: eor ! end of record indicator
    integer(kind=ikind) :: total_cols = 0 ! initialize total number of columns 
    integer(kind=ikind) :: total_columns_in_file ! variable that would be returned
    integer(kind=ikind) :: setflag = 0 ! to know whether it is space delimited or comma delimited
    !***************************************************************************
    total_cols = 0 ! initialize total columns to zero
    setflag = 0 ! initilize setflag to zero
    eor = 0 ! initialize end of record indicator
    name=trim(adjustl(inp_direc))//trim(adjustl(file_name))!file name with directory
    ! open the specified file
    open(unit=1,file = name,form="formatted",status="old",action="read")
    do while (eor == 0)
        read(1,"(a1)",advance="no", iostat=eor) buffer
        setflag = setflag+1;
        if (setflag==1 .and. (ichar(buffer)==ichar(",") .or. ichar(buffer)==ichar(" ") & 
        .or. ichar(buffer)== ichar("	"))) then
            read(1,"(a1)",advance="no", iostat=eor) buffer
        endif
        if (ichar(buffer)==ichar(".") .or. (ichar(buffer)>=ichar("0")& 
        .and. ichar(buffer)<=ichar("9"))) then
            do while (ichar(buffer)==ichar(".") .or. (ichar(buffer)>=ichar("0") & 
            .and. ichar(buffer)<=ichar("9")))
                read(1,"(a1)",advance="no", iostat=eor) buffer
            enddo
            total_cols=total_cols+1;
        elseif (ichar(buffer)<ichar("0").or. ichar(buffer)>ichar("9")) then
            do  
                if (ichar(buffer)==ichar(",") .or. ichar(buffer)==ichar(" ") & 
                .or. ichar(buffer)== ichar("	")) then
                    exit
                endif
                read(1,"(a1)",advance="no", iostat=eor) buffer
            enddo
            total_cols=total_cols+1
        endif
    enddo
    close(1)
    total_columns_in_file = total_cols
    end function total_columns

    function file_write_array(out_direc,file_name,array,rows,cols) result (status)
    !**************************************************************************
    !written by : vineet yadav  									   
    !date : 8/29/2010
    !purpose : write space delimited text file
    !***********specification part*********************************************
    implicit none
    character (len=*), intent(in):: out_direc
    character (len=*), intent(in):: file_name
    character (len=200):: name ! combined name of out_direc and file_name
    integer(kind=ikind) :: rows !number of lines in the file
    integer(kind=ikind) :: cols !number of cols in the file
    integer(kind=ikind) :: i ! loop iterator
    logical :: status ! returns true if successful 
    real(kind=rkind), dimension(:,:) :: array(rows,cols)
    !**********execution part**************************************************
    ! open the specified file
    name=trim(adjustl(out_direc))//trim(adjustl(file_name))
    open(unit=1,file=name,form="formatted",status="replace",action="write")
    ! write_lines
    do i=1,rows
        write(1,*) array(i,:)
    enddo
    close(1)
    status = .true.
    end function file_write_array

    subroutine error_check(iocheck,fname)
    !**************************************************************************
    !written by : Mae Qiu
    !date : 7/12/13
    !purpose : give error messages if the read function fails
    !***********specification part*********************************************
    implicit none
    character (len=*) :: fname
    integer(kind=ikind),intent(in):: iocheck
    if (iocheck <0) then
       write(*,*) 'End of File or End of Record reached in ',fname
       write(*,*) '(you tried to read more numbers than existed in the record)'
       stop
    elseif (iocheck >0) then
       write(*,*) 'Error opening or unrecoverable error when reading ',fname
       stop
    endif
    
    end subroutine error_check
    end module file_io
