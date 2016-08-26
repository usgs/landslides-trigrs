c !  subroutine to read an ascii grid file and store in a 1-d array
c !  by Rex L. Baum, USGS May 2001 latest revison 13 Mar 2013, RLB
c !  single precision
c !
	subroutine srdgrd(grd,pth,ncol,nrow,celsiz,nodat,
     +	pf,pf1,ctr,imax,temp,u,infil,param,header,u1)
	implicit none
	integer grd,pth,i,m,ncol,nrow,ctr,imax,u,u1,lncnt
	double precision param(6),nodat,celsiz,cns,cew
	double precision east,west,north,south
	real pf(imax),pf1(grd),temp(pth),nodats
	character*14 header(6),h1
	character*255 infil
c !  	
	open(u,file=trim(infil),status='old',err=23)
	do 200, m=1,6
        lncnt=m
	read(u,*,err=130) header(m),param(m)
	header(m)=adjustl(header(m))
  200	continue
c ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
  	nodat=-9999.d0
  	celsiz=-10.d0
	do 210, m=1,6	   
	h1=header(m); h1=adjustl(h1)
	if(trim(h1).eq.'ncols'.or.trim(h1).eq.'NCOLS')ncol=int(param(m))
	if(trim(h1).eq.'nrows'.or.trim(h1).eq.'NROWS')nrow=int(param(m))
	if(trim(h1).eq.'cellsize'.or.trim(h1).eq.'CELLSIZE')
     + celsiz=param(m)
	if(trim(h1).eq.'NODATA_value')nodat=param(m)
	if(trim(h1).eq.'nodata_value'.or.trim(h1).eq.'NODATA_VALUE') 
     + nodat=param(m)
	if(trim(h1).eq.'cols:'.or.trim(h1).eq.'COLS:')ncol=int(param(m))
	if(trim(h1).eq.'rows:'.or.trim(h1).eq.'ROWS:')nrow=int(param(m))
	if(trim(h1).eq.'east:'.or.trim(h1).eq.'EAST:')east=param(m)
	if(trim(h1).eq.'west:'.or.trim(h1).eq.'WEST:')west=param(m)
	if(trim(h1).eq.'north:'.or.trim(h1).eq.'NORTH:')north=param(m)
	if(trim(h1).eq.'south:'.or.trim(h1).eq.'SOUTH:')south=param(m)
  210	continue
  	if (celsiz.le.0) then
  	  cew=abs(east-west)/ncol
  	  cns=abs(north-south)/nrow
  	  if (cew.eq.cns) then
  	    celsiz=cew
  	  else
  	    celsiz=sqrt(cew*cns)
  	    write(*,*) 'Rectangular cells ',cew, ' X ', cns
  	    write(u1,*) 'Rectangular cells ',cew, ' X ', cns
  	  end if
  	end if
  	if (ncol*nrow .gt. grd) then
  	  write(*,*)'Grid file exceeds array size in subroutine srdgrd'
	  write (*,*)'--> ',trim(infil)
  	  write(*,*)'Check intialization file row and column values.'
  	  write(u1,*)'Grid file exceeds array size in subroutine srdgrd'
	  write (u1,*)'--> ',trim(infil)
  	  write(u1,*)'Check intialization file row and column values.'
      	  write(*,*)'Press RETURN to exit'
  	  read*
	  close(u)
	  close(u1)
	  stop
  	end if
  	nodats=nodat
	ctr=0
	do 120, m=1,nrow
c !  next sequence of lines read data in but skips no_data values
c !  count maintained by ctr should coincide with node numbers from GIS
c !  pf1() keeps track of positions of nodata values so that results
c !  can be written out in grid format.
        lncnt=m+6
	read(u,*,end=125,err=130) (temp(i), i=1,ncol)
	do 250, i=1,ncol
	  pf1(i+(m-1)*ncol)=temp(i)
	  if(temp(i).ne.nodats) then
	    ctr=ctr+1
	    if (ctr>imax) then
              write (*,* )'Subroutine srdgrd reports'
  	      write(*,*) 'Number of data cells exceeds array size'
	      write (*,*) '--> ',trim(infil)
  	      write(*,*) 'Check imax value in grid_size.txt file and
     + no-data values in input grid.'	    	
              write (u1,*) 'Subroutine srdgrd reports'
  	      write(u1,*) 'Number of data cells exceeds array size'
	      write (u1,*) '--> ',trim(infil)
  	      write(u1,*) 'Check imax value in grid_size.txt file and
     + no-data values in input grid.'	    	
  	      write(*,*) 'Press RETURN to exit'
  	      read*
  	      close(u)
	      close(u1)
	      stop 'Error at line 77 in subroutine srdgrd'
	    end if
	    pf(ctr)=temp(i)
	  end if
  250	continue
  120	continue
  125	close(u)
	return
   23	continue
   	write(*,*)'***Error opening input file in subroutine srdgrd***'
	write(*,*)'--> ',trim(infil)
	write(*,*)'Check file name and location'
   	write(u1,*)'***Error opening input file in subroutine srdgrd***'
	write(u1,*)'--> ',trim(infil)
	write(u1,*)'Check file name and location'
  	write(*,*)'Press RETURN to exit'
  	read*
	close(u)
	close(u1)
	stop '-23 in srdgrd()'
 130	continue
  	write (*,*) 'Error reading grid file, line '
	write(*,*)'--> ',trim(infil), lncnt
  	write (u1,*) 'Error reading grid file, line '
	write(u1,*)'--> ',trim(infil), lncnt
        write(*,*) 'Press RETURN to exit'
        read*
	close(u)
	close(u1)
	stop '-130 in subroutine srdgrd'
	end

