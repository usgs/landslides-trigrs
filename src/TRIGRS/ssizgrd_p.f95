!  subroutine to read an ascii grid file (elevations) and determine its size (rows, columns, & data cells)
!  by Rex L. Baum, USGS Feb 2011 latest revison 13 Mar 2013
!  single precision
!
! MPI calls & MPI-related changes by M. Alvioli, November 2014  
!
subroutine ssizgrd_p(row,col,celsiz,nodat,ctr,u,infil,header,ulog)
  use mpi
  implicit none
  integer,parameter:: double=kind(1d0)
  integer ::i,m,ctr,u,ulog 
  integer :: col,row,ctall
  real (double):: nodat,celsiz,cns,cew,param(6)
  real (double):: east,west,north,south
  real:: nodats
  real, allocatable:: temp(:)
  character (len=14):: header(6)
  character (len=255):: infil
  integer myrank,isize,ierr
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  !  	
  infil=adjustl(infil) ! adjustl statements added to improve compatibility with other compilers 14 Feb 2013 RLB
  if(myrank.eq.0) then
     open(u,file=trim(infil),status='old',err=23)
     do m=1,6
        read(u,*) header(m),param(m)
        header(m)=adjustl(header(m))
     end do
  endif
  do m=1,6
     Call MPI_BCAST(header(m),14,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
  enddo
  Call MPI_BCAST(param,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
  nodat=-9999.d0
  celsiz=-10.d0
  do m=1,6
     if (trim(header(m)).eq.'ncols') col=int(param(m))
     if (trim(header(m)).eq.'nrows') row=int(param(m))
     if (trim(header(m)).eq.'cellsize') celsiz=param(m)
     if (trim(header(m)).eq.'NODATA_value') nodat=param(m)
     if (trim(header(m)).eq.'nodata_value') nodat=param(m)
     if (trim(header(m)).eq.'cols:') col=int(param(m))
     if (trim(header(m)).eq.'rows:') row=int(param(m))
     if (trim(header(m)).eq.'east:') east=param(m)
     if (trim(header(m)).eq.'west:') west=param(m)
     if (trim(header(m)).eq.'north:') north=param(m)
     if (trim(header(m)).eq.'south:') south=param(m)
  end do
  if (celsiz.le.0) then
     cew=abs(east-west)/col
     cns=abs(north-south)/row
     if (cew.eq.cns) then
        celsiz=cew
     else
        celsiz=sqrt(cew*cns)
        if(myrank.eq.0) then
           write(*,*) 'Rectangular cells ',cew, ' X ', cns
           write(ulog,*) 'Rectangular cells ',cew, ' X ', cns
        endif
     end if
  end if
  nodats=nodat
  allocate(temp(col))
  ctr=0; ctall=0
  row_loop: do m=1,row
     !  next sequence of lines read data in but skips no_data values
     !  count maintained by ctr should coincide with node numbers from GIS
     if(myrank.eq.0) read(u,*,end=125) (temp(i), i=1,col) 
     Call MPI_BCAST(temp,col,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
     col_loop: do i=1,col
        ctall=ctall+1
        if(temp(i) /= nodats) then
           ctr=ctr+1
        end if
     end do col_loop
  end do row_loop
125 close(u)
  if(myrank.eq.0) then
     write(*,*) ctr,' = number of data cells'
     write(*,*) ctall,' = total number of cells'
  endif
  return
23 continue
  if(myrank.eq.0) then
     write (*,*) '*** Error opening input file in subroutine ssizgrd ***'
     write (*,*) '--> ',trim(infil)
     write (*,*) 'Check file name and location'
     write (ulog,*) '*** Error opening input file in subroutine ssizgrd ***'
     write (ulog,*) '--> ',trim(infil)
     write (ulog,*) 'Check file name and location'
     close(u)
     close(ulog)
  endif
  call MPI_FINALIZE(ierr)
  stop '- Error in ssizgrd()'
end subroutine ssizgrd_p
	
