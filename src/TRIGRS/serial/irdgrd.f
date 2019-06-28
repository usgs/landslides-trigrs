c !  subroutine to read an ascii grid file of integer values
c !  and store in a 1-d array
c ! Rex L. Baum, USGS, spring 2001
c !  Dec 2003, added code to handle GRASS GIS ASCII grid files
c ! 
	subroutine irdgrd(grd,pth,ncol,nrow,celsiz,inodat,mnd,
     +	y,y1,ctr,imax,itemp,u,infil,param,header,u1)
	implicit none
	integer grd,pth,i,m,ncol,nrow,ctr,imax,u,inodat,u1
	integer y(imax),y1(grd),itemp(pth),mnd,lncnt
	double precision param(6),celsiz,cns,cew,nodat
	double precision east,west,north,south
	character*14 header(6),h1
	character*255 infil
c !  	
        infil=adjustl(infil)
	open(u,file=trim(infil),status='old',err=10)
	do 200, m=1,6
	lncnt=m
	read(u,*,err=130) header(m),param(m)
	header(m)=adjustl(header(m))
  200	continue
c ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
  	inodat=-9999
  	celsiz=-10.
	do 210, m=1,6
	h1=header(m); h1=adjustl(h1)
	if(trim(h1).eq.'ncols'.or.trim(h1).eq.'NCOLS')ncol=int(param(m))
	if(trim(h1).eq.'nrows'.or.trim(h1).eq.'NROWS')nrow=int(param(m))
        if(trim(h1).eq.'cellsize'.or.trim(h1).eq.'CELLSIZE')
     + celsiz=param(m)
	if(trim(h1).eq.'cols:'.or.trim(h1).eq.'COLS:')ncol=int(param(m))
	if(trim(h1).eq.'rows:'.or.trim(h1).eq.'ROWS:')nrow=int(param(m))
	if(trim(h1).eq.'east:'.or.trim(h1).eq.'EAST:')east=param(m)
	if(trim(h1).eq.'west:'.or.trim(h1).eq.'WEST:')west=param(m)
	if(trim(h1).eq.'north:'.or.trim(h1).eq.'NORTH:')north=param(m)
	if(trim(h1).eq.'south:'.or.trim(h1).eq.'SOUTH:')south=param(m)
	if(trim(h1).eq.'NODATA_value')  then
	  nodat=(param(m))
	  mnd=m
	  inodat=int(nodat)
	end if
	if(trim(h1).eq.'nodata_value'.or.trim(h1).eq.'NODATA_VALUE')then
	  nodat=(param(m))
	  mnd=m 
	  inodat=int(nodat)
	end if
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
  	  write(*,*)'Grid file exceeds array size in subroutine irdgrd'
	  write (*,*)'--> ',trim(infil)
  	  write(*,*)'Check intialization file row and column values.'
  	  write(u1,*)'Grid file exceeds array size in subroutine irdgrd'
	  write (u1,*)'--> ',trim(infil)
  	  write(u1,*)'Check intialization file row and column values.'
  	  write(*,*)'Press RETURN to exit'
  	  read*
	  close(u)
	  close(u1)
	  stop '-11 in subroutine irdgrd'
  	end if
	ctr=0
	do 120, m=1,nrow
c !  next sequence of lines read data in but skips no_data values
c !  count maintained by ctr should coincide with node numbers from GIS
c !  y1() keeps track of positions of nodata values so that results
c !  can be written out in grid format.
        lncnt=m+6
  	read(u,*,end=125,err=130) (itemp(i), i=1,ncol)
	do 250, i=1,ncol
	  y1(i+(m-1)*ncol)=itemp(i)
	  if(itemp(i).ne.inodat) then
	  ctr=ctr+1
	  if (ctr>imax) then
            write (*,* )'Subroutine irdgrd reports'
  	    write(*,*) 'Number of data cells exceeds array size'
	    write (*,*) '--> ',trim(infil)
  	    write(*,*) 'Check imax value in grid_size.txt file and 
     + no-data values in input grid.'
            write (*,* )'Subroutine irdgrd reports'
  	    write(u1,*) 'Number of data cells exceeds array size'
	    write (u1,*) '--> ',trim(infil)
  	    write(u1,*) 'Check imax value in grid_size.txt file 
     + no-data values in input grid.'
  	    write(*,*) 'Press RETURN to exit'
  	    read*
  	    close(u)
	    close(u1)
	    stop 'Error in subroutine irdgrd'
	  end if
	  y(ctr)=itemp(i)
	  end if
  250	continue
  120	continue
  125	close(u)
	return
   10	continue
   	write(*,*)'***Error opening input file in subroutine irdgrd***'
	write(*,*)'--> ',trim(infil)
	write(*,*)'Check file name and location'
   	write(u1,*)'***Error opening input file in subroutine irdgrd***'
	write(u1,*)'--> ',trim(infil)
	write(u1,*)'Check file name and location'
  	write(*,*)'Press RETURN to exit'
  	read*
	close(u)
	close(u1)
	stop '-10 in subroutine irdgrd'
  130	continue
  	write (*,*) 'Error reading grid file, line '
	write(*,*)'--> ',trim(infil), lncnt
  	write (u1,*) 'Error reading grid file, line '
	write(u1,*)'--> ',trim(infil), lncnt
                write(*,*) 'Press RETURN to exit'
                read*
	close(u)
	close(u1)
	stop '-130 in subroutine irdgrd'
	end
