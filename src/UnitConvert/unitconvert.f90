!  Program to compute unit conversions of grids
!  
!  Rex L. Baum, USGS, spring 2001
!   conversion to Fortran 90, March 2002
	program unitconvert
! 
	implicit none
 	integer:: grd,row,col
	integer:: i 
	integer,parameter:: double=kind(1d0)
  	integer:: sctr 
	integer:: ncol,nrow,u0,u1,u2,u3
	logical:: ans
	real, allocatable:: x(:),pf1(:),temp(:)
	real(double):: nodat,fac 
	real(double):: celsiz,param(6),ti
	real:: nodata
	character:: tb*1 
  	character (len=255):: infil,outfil,heading,title
	character (len=14):: header(6)
	character (len=9):: init
	character (len=8):: date
	character (len=10):: time
	character (len=6):: vrsn
	character (len=11):: bldate
	data u0,u1,u2,u3/10,11,12,13/
 	tb=char(9)
     	call date_and_time(date,time)
	vrsn='1.0.03'; bldate='24 Mar 2016'
! Display program information	
	write (*,*) '   UnitConvert: A simple utility for'
	write (*,*) '  converting ASCII grid data from one'
	write (*,*) '      system of units to another'
	write (*,*) '       Version ', vrsn,', ',bldate
	write (*,*) '             By Rex L. Baum'
	write (*,*) '        U.S. Geological Survey'
	write (*,*) '-----------------------------------------'
	write (*,*) ''
!  Open log file
	outfil='UnitConvertLog.txt'; outfil=adjustl(outfil)
	open (u1,file=trim(outfil),status='unknown',err=20)
	write (u1,*) ''
     	write (u1,*) 'Starting UnitConvert'
     	write (u1,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
	write (u1,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
!  Open initialization file
	init='uc_in.txt' ; init=adjustl(init)
	inquire (file=init,exist=ans)
	if(ans)then
 	  open(u3,file=init,status='old',err=10)
	else
	  write (*,*) 'Cannot locate default initialization file'
	  write (*,*) 'Type name of initialization file and'
	  write (*,*) 'press RETURN to continue'
	  read (*,'(a)') init; init=adjustl(init)
	  open (u3,file=trim(init),status='old',err=10)
	end if
!  Read initialization file and echo in log file
	write (u1,*) 'Initialization file -->', trim(init)
	write (u1,*) '-- LISTING OF INITIALIZATION FILE --'	
	read (u3,'(a)',err=15) heading; heading=adjustl(heading)
	read (u3,'(a)',err=15) title; title=adjustl(title)
  	write (u1,*) trim(heading)
  	write (u1,*) trim(title)
	read (u3,'(a)',err=15) heading; heading=adjustl(heading)
	read (u3,*,err=15) row, col
  	write (u1,*) trim(heading)
  	write (u1,*) row, col
	read (u3,'(a)',err=15) heading; heading=adjustl(heading)
	read (u3,'(a)',err=15) infil; infil=adjustl(infil)
  	write (u1,*) trim(heading)
  	write (u1,*) trim(infil)
	read (u3,'(a)',err=15) heading; heading=adjustl(heading)
	read (u3,'(a)',err=15) outfil; outfil=adjustl(outfil)
  	write (u1,*) trim(heading)
  	write (u1,*) trim(outfil)
	read (u3,'(a)',err=15) heading; heading=adjustl(heading)
	read (u3,*,err=15) fac
  	write (u1,*) trim(heading)
  	write (u1,*) fac
  	write (u1,*) '-- END INITIALIZATION DATA --'
  	write (u1,*) ''
	close(u3)
	grd=row*col
	allocate (x(grd),pf1(grd),temp(col))
!  Initialize arrays
	x=0.
	pf1=0.
!  Read gridded data from GIS
	call srdgrd(grd,col,ncol,nrow,celsiz,nodat, &
        & x,pf1,sctr,grd,temp,u0,infil,param,header,u1)
	nodata=nodat
     	do i=1,sctr
	  x(i)=x(i)*fac
	end do
!  Write output grid
	ti=tiny(param(1))
	call ssvgrd(x,grd,pf1,nrow,ncol,u2,nodata,param,u1, &
	& outfil,ti,header)
	write (u1,*) 'UnitConvert finished normally'
     	call date_and_time(date,time)
     	write (u1,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
	write (u1,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
     	close (u1)
   	stop
! Error handling if unable to open or read files	
   10	continue
	write (*,*) '*** Error opening initialization file ***'
	write (*,*) '--> ',trim(init)
	write (*,*) 'Check file name and location'
  	write (u1,*) 'Error reading initialization file'
	write (u1,*) '--> ',trim(init)
	write (u1,*) 'Check file contents and organization'
	pause 'Press RETURN to exit'
	stop
   15	continue
  	write (*,*) 'Error reading initialization file'
	write (*,*) '--> ',trim(init)
	write (*,*) 'Check file contents and organization'
  	write (u1,*) 'Error reading initialization file'
	write (u1,*) '--> ',trim(init)
	write (u1,*) 'Check file contents and organization'
	pause 'Press RETURN to exit'
	stop
   20	continue
  	write (*,*) 'Error opening output file'
	write (*,*) '--> ',trim(outfil)
	write (*,*) 'Check file path and status'
  	write (u1,*) 'Error opening output file'
	write (u1,*) '--> ',trim(outfil)
	write (u1,*) 'Check file path and status'
	pause 'Press RETURN to exit'
	stop
   	end program unitconvert
