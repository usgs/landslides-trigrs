!  Program to check nodata cell locations in ascii grids
!  Log file contains results
!  Rex L. Baum, USGS, Fall 2002
	program gridmatch
! 
	implicit none
 	integer:: grd,row,col,patlen
	integer:: i,j,num,curcol,currow
	integer,parameter:: DOUBLE=kind(1D0)
  	integer:: sctr,sctr1,mmctr,mmctr1,mmctr2,clcnt,dsctr
	integer:: ncol,nrow,ncol1,nrow1,u(6),cel_ctr
	logical:: ans,lskip
	real, allocatable:: mx(:),x(:),x1(:),pf1(:),pf2(:),temp(:)
	real(DOUBLE):: nodat,ti
	real(DOUBLE):: celsiz,param(6)
	real:: nodata,nodata1 
	character:: tb*1
  	character (len=255),allocatable:: infil(:)
  	character (len=255):: outfil,init,heading,demfil,scratch
  	character (len=224):: outpath,outpath1
  	character (len=255)::fname
	character (len=14):: header(6)
	character (len=10):: stime
	character (len=11):: bldate
	character (len=8):: sdate
	character (len=6):: vrsn
     	call date_and_time(sdate,stime)
	vrsn='1.0.05'; bldate='24 Mar 2016'
	u=(/11,12,13,14,15,16/)
 	tb=char(9)
!  Open log file
	outfil='GridMatchLog.txt'; outfil=adjustl(outfil)
	open (u(2),file=trim(outfil),status='unknown',err=20)
	write (u(2),*) ''
     	write (u(2),*) 'Starting GridMatch'
	write (u(2),*) 'Version ', vrsn,' ',bldate
     	write (u(2),*) 'Date ',sdate(5:6),'/',sdate(7:8),'/',sdate(1:4)
	write (u(2),*) 'Time ',stime(1:2),':',stime(3:4),':',stime(5:6)
! Display program information	
	write (*,*) ''
	write (*,*) '   GridMatch: A simple utility for'
	write (*,*) '    checking ASCII grid data and'
	write (*,*) '      comparing the locations of'
	write (*,*) '             no-data  cells'
	write (*,*) '       Version ', vrsn,' ',bldate
	write (*,*) '             By Rex L. Baum'
	write (*,*) '        U.S. Geological Survey'
	write (*,*) '-----------------------------------------'
	write (*,*) ''
!  Open initialization file
	init='gm_in.txt' ; init=adjustl(init)
	inquire (file=trim(init),exist=ans)
	if(ans)then
 	  open(u(4),file=trim(init),status='old',err=10)
	else
	  write (*,*) 'Cannot locate default initialization file'
	  write (*,*) 'Type name of initialization file and'
	  write (*,*) 'press RETURN to continue'
	  read (*,'(a)') init
	  open (u(4),file=trim(init),status='old',err=10)
	end if
!  Read initialization file and echo in log file
	write (u(2),*) 'Initialization file -->', trim(init)
	write (u(2),*) '-- LISTING OF INITIALIZATION FILE --'	
	read (u(4),'(a)',err=15) heading; heading=adjustl(heading)
	read (u(4),*,err=15) num 
  	write (u(2),*) trim(heading)
  	write (u(2),*) num
	allocate (infil(num))
	  read (u(4),'(a)',err=15) heading; heading=adjustl(heading)
  	  write (u(2),*) trim(heading)
	do i=1,num
	  read (u(4),'(a)',err=15) infil(i)
	  infil(i)=adjustl(infil(i))
  	  write (u(2),*) trim(infil(i))
	end do
        demfil=infil(1)
        patlen=scan(demfil,'/\',.true.)
        outpath=demfil(1:patlen); outpath=adjustl(outpath) ! path to elevation grid
	write(u(2),*) 'Path to elevation grid and output files'
	write(u(2),*) trim(outpath)
  	write (u(2),*) '-- END INITIALIZATION DATA --'
  	write (u(2),*) ''
	close(u(4))
! Obtain number of rows, columns and data cells from the master grid
	call ssizgrd(row,col,celsiz,nodat,clcnt,u(1),infil(1),header,u(2))
	grd=row*col
	allocate (mx(grd),x(grd),x1(grd),pf1(grd),temp(col))
	allocate(pf2(grd))
!  Initialize arrays
	mx=0.
	pf1=0.
	x=0.;x1=0
	pf2=0.
!  Read master grid data from GIS
	call srdgrd(grd,col,ncol,nrow,celsiz,nodat, &
        & mx,pf1,sctr,grd,temp,u(1),infil(1),param,header,u(2))
	nodata=nodat
	write (u(2),*) 'Reference grid,'
	write (u(2),*) trim(infil(1))
	write (u(2),*) 'contains ',sctr,' data cells'
	write (u(2),*) 'contains ',ncol,' columns'
	write (u(2),*) 'contains ',nrow,' rows'
	write (u(2),*) '--------------------------------------------'
!  Read check, and adjust other grid data from GIS
     	do i=2,num
     	  lskip=.false.
	  call srdgrd(grd,col,ncol1,nrow1,celsiz,nodat, &
          & x,pf2,sctr1,grd,temp,u(3),infil(i),param,header,u(2))
	  nodata1=nodat
	  write(*,*) 'Testing ', trim(infil(i))
	  write (u(2),*) ''
	  write (u(2),*) 'Results for grid,'
	  write (u(2),*) trim(infil(i))
	  write (u(2),*) 'contains ',sctr1,' data cells'
	  if(sctr1/=sctr) then
	    write (u(2),*) '**** Does not match test grid! ****'
	  end if
	  write (u(2),*) 'contains ',ncol1,' columns'
	  if(ncol1/=ncol) then
	    write (u(2),*) '**** Does not match test grid! ****'
	    lskip=.true.
	  end if
	  write (u(2),*) 'contains ',nrow1,' rows'
	  if(nrow1/=nrow) then
	    write (u(2),*) '**** Does not match test grid! ****'
	    lskip=.true.
	  end if
	  if(lskip) then
	    write(u(2),*) 'Grid shape or size mismatch, skipping &
	     &comparison of nodata cells'
	    cycle
	  end if
	  mmctr1=0; mmctr2=0 ! Scan grids to find mismatching cells
	  do j=1,grd
	    if (pf1(j)==nodata) then
	      if (pf2(j)/=nodata1) then
	        curcol=mod(j,ncol)
	        currow=1+j/ncol
	        if(curcol==0) then
	          curcol=ncol
	          currow=currow-1
	        end if
	        write(u(2),*) 'Cell mismatch found! Cell',j
	        write(u(2),*) 'at row, ',currow, 'column, ',curcol
	        write(u(2),*) 'is a nodata cell in the reference grid'
	        write(u(2),*) 'and a data cell in grid ', trim(infil(i))
	        mmctr1=mmctr1+1
	      end if 
	     end if
	      if (pf1(j)/=nodata) then
	        if (pf2(j)==nodata1) then
	          write(u(2),*) 'Cell mismatch found! Cell',j
	          write(u(2),*) 'at row, ',1+j/ncol, 'column, ',mod(j,ncol)
	          write(u(2),*) 'is a data cell in the reference grid'
	          write(u(2),*) 'and a nodata cell in grid ', trim(infil(i))
	          mmctr2=mmctr2+1
	        end if
	      end if
	  end do
	  mmctr=mmctr1+mmctr2
	  write(*,*) 'Number of mismatches found: ',mmctr1, ' + ', mmctr2, ' = ', mmctr
	  write(u(2),*) 'Number of mismatches found: ',mmctr1, ' + ', mmctr2, ' = ', mmctr
	  write (u(2),*) '--------------------------------------------'
	if(mmctr1>0 .and. mmctr2==0) then
	  write(*,*) 'Saving corrected grid with matching no-data cells'
	  cel_ctr=0
                   do j=1,grd
                      if(pf1(j)==nodata) cycle
                      cel_ctr=cel_ctr+1
                      x1(cel_ctr)=x(j)
                   end do
! name corrected output grid
                  write(scratch,*) infil(i);scratch=trim(adjustl(scratch))
                  patlen=scan(scratch,'/\',.true.)
                  outpath1=scratch(1:patlen); outpath1=adjustl(outpath1) ! path to ith input grid
                  fname=scratch(patlen+1:255)
!                  write(*,*) trim(scratch),patlen,trim(outpath1),trim(fname)
	  outfil=trim(outpath1)//'GM'//trim(fname)
	  outfil=adjustl(outfil)
	  write(*,*) 'Saving ', trim(outfil)
	  ti=tiny(param(1))
                 call ssvgrd(x1,grd,pf1,nrow,ncol,u(6),nodata1,param,u(2),outfil,ti,header)
                 end if
	  write (*,*) '--------------------------------------------'
	end do
! write grid size parameters for use by TRIGRS RB, 3/18/11 
	if(mmctr == 0) then
	  outfil=trim(outpath)//'GMgrid_size.txt'
	  outfil=adjustl(outfil)
          open (u(5),file=trim(outfil),status='unknown',err=20)
	  write (u(5),*) 'imax      nrow      ncol      nwf'
	  dsctr=1 ! dsctr is computed by TopoIndex; dsctr=1 is default value for no runoff routing.
	  write (u(5),*) sctr,nrow,ncol,dsctr
     	  write (u(5),*) ''
     	  close (u(5))
	end if
! RB
	write (u(2),*) 'GridMatch finished normally'
     	call date_and_time(sdate,stime)
     	write (u(2),*) 'Date ',sdate(5:6),'/',sdate(7:8),'/',sdate(1:4)
	write (u(2),*) 'Time ',stime(1:2),':',stime(3:4),':',stime(5:6)
     	close (u(2))
   	stop
! Error handling if unable to open or read files	
   10	continue
	write (*,*) '*** Error opening initialization file ***'
	write (*,*) '--> ',trim(init)
	write (*,*) 'Check file name and location'
  	write (u(2),*) 'Error reading initialization file'
	write (u(2),*) '--> ',trim(init)
	write (u(2),*) 'Check file contents and organization'
	pause 'Press RETURN to exit'
	stop
   15	continue
  	write (*,*) 'Error reading initialization file'
	write (*,*) '--> ',trim(init)
	write (*,*) 'Check file contents and organization'
  	write (u(2),*) 'Error reading initialization file'
	write (u(2),*) '--> ',trim(init)
	write (u(2),*) 'Check file contents and organization'
	pause 'Press RETURN to exit'
	stop
   20	continue
  	write (*,*) 'Error opening output file'
	write (*,*) '--> ',trim(outfil)
	write (*,*) 'Check file path and status'
  	write (u(2),*) 'Error opening output file'
	write (u(2),*) '--> ',trim(outfil)
	write (u(2),*) 'Check file path and status'
	pause 'Press RETURN to exit'
	stop
   	end program gridmatch
