!  Routine to create file headers for various list file formats (1/20/2011)
!
!  Rex L. Baum, 21 Apr 2010, latest revision 07 Feb 2012
	subroutine prpijz(u1,ulog,profil,ncol,nrow,header,vrsn)
	use grids; use input_vars 
	use model_vars
	use input_file_defs
	implicit none
	integer::i,j,jj,ctr,ncol,nrow,ulog,u1,ctr1,swt,m,nvar,outrow,zval
	real (double):: test0
	real (double):: xmin,xmax,ymin,ymax,zlo,zhi,plo,pmax,th_min,th_max
	character (len=255)::outfil
	character (len=18):: profil
	character (len=14):: header(6)
	character (len=8)::fini(2)
	character (len=7):: vrsn ! added 12/28/2010
	character (len=4)::stp
	fini=(/'Infinite','Finite  '/)
	swt=1
	if(mmax>0) swt=2
	fini(swt)=adjustl(fini(swt))
! Open file and write header
! list file for storing depth profile on each cell
	if (flag<=-1 .and. flag>=-3) then
	  outfil=trim(folder)//trim(profil)//trim(suffix)//'.txt'
	  outfil=adjustl(outfil)
	  open (u1,file=trim(outfil),status='unknown',err=10)
	  write(u1,'(a,a)',err=11) 'Depth profiles at each cell, TRIGRS, v. ',vrsn
	  write(u1,'(a,a)',err=11) trim(fini(swt)),'-depth no-flow boundary' ! format corrected 5/3/10 RLB
	  write(u1,'(a)',err=11) 'Cell Number, Slope angle, Step#, Time'
 	  if (flag==-1) write(u1,'(a)',err=11) 'Z         P         FS'
  	  if (flag==-2) write(u1,'(a)',err=11) &
          &'Z         P      Pzero     Ptran      Pbeta       FS'
  	  if (flag==-3 .and. unsat0) then
  	    write(u1,'(a)',err=11) 'Z         P         FS      TH' !added 4/14/2010 RLB
  	  else if (flag==-3) then
  	    write(u1,'(a)',err=11) 'Z         P         FS'
  	  end if
  	  return ! Reorganized if block and added return statement 9/6/2011, RLB
	end if
! ijz file for depth profile of each cell
  	if(flag==-4 .or. flag==-5 .or. flag==-6) then ! revised 12/23/2010
	    if(deepz>0) then
	      zval=nzs+3
	    else
	      zval=nzs+2
	    end if
  	  profil='TR_ijz_p_th_'
  	  profil=adjustl(profil)
	  do j=1,nout ! added loop to create series of files 15 Nov 2010 RLB
	    uijz(j)=j-1+uijz(1) !12/6/2010, added consecutive list of file unit numbers.
	    write(stp,'(i4)') j
	    stp=adjustl(stp)
	    outfil=trim(folder)//trim(profil)//trim(suffix)//'_'//trim(stp)//'.txt'
            outfil=adjustl(outfil)
	    open(uijz(j),file=trim(outfil),status='unknown',err=10) ! Revised formating in succeeding write statements, RLB 9/2/2011
	    write(uijz(j),'(a,a)',err=11) '# Comments: Quasi-3-D pressure head data from TRIGRS, v. ',vrsn
	    write(uijz(j),'(a,a)',err=11) '# ', title ! revised format 1/6/2011 RLB
	    if (unsat0) then
	      write(uijz(j),'(a)',err=11) '# i  j  z  p  th'
	    else
	      write(uijz(j),'(a)',err=11) '# i  j  z  p'
	    end if
	    write(uijz(j),'(a12,1x,i12)',err=11) '# timestep= ', ksav(j)
	    write(uijz(j),'(a8,1x,g20.7)',err=11) '# time= ', tsav(j)
! New SCOOPS header format, updated 7/27/2011 RLB, Revised 2/28/2012 RLB
	    write(uijz(j),'(a)',err=11) 'coords'
	    write(uijz(j),'(a)',err=11) 'ijz'
	  end do
	end if
	ctr=0;ctr1=0
	if(flag<=-4) then
!  Map i-j values to cell numbers for flag=-4, -5, -6
	  do j=1,nrow
	    do i=1,ncol
	      ctr1=ctr1+1
	      test0=abs(pf1(ctr1)-test1) ! check against default no-data value
  	      if (test0.le.0.1) then
  	        cycle
  	      else
  	        ctr=ctr+1
	        jj=nrow-j+1 ! transform from origin at upper left to lower left corner of grid
	        ix(ctr)=i; jy(ctr)=jj
  	      end if
	    end do
	  end do
  	end if
! Copy coordinates of SW corner of grid for xyz file formats
  	if (flag<=-7 .and. flag>=-9) then
	  do m=1,6
	    if (trim(header(m)).eq.'xllcorner') xllc=param(m)
	    if (trim(header(m)).eq.'yllcorner') yllc=(param(m))
	    if (trim(header(m)).eq.'west:') xllc=param(m)
	    if (trim(header(m)).eq.'south:') yllc=param(m)
	  end do
  	end if
!
	if(flag<=-7 .and. flag>=-9) then ! Added 7/27/2011, RLB, Xmdv format, readable by VisIt 
            profil='TR_xyz_p_th_'
            profil=adjustl(profil)
            nvar=5
! xmin, ymin at lower left corner of lower left cell, xmax, ymax at upper right corner of upper right grid cell
  	  xmin=xllc
  	  xmax=xllc+dfloat(ncol)*celsiz
  	  ymin=yllc
  	  ymax=yllc+dfloat(nrow)*celsiz
  	  pmax=maxval(zmax)
  	  plo=-pmax
  	  if(deepz > pmax) pmax=deepz
  	  zlo=minval(elev)-pmax
  	  zhi=maxval(elev)
  	  th_min=0.
  	  th_max=maxval(ths)
	  do j=1,nout
	    uijz(j)=j-1+uijz(1)
	    write(stp,'(i4)') j
	    stp=adjustl(stp)
	    outfil=trim(folder)//trim(profil)//trim(suffix)//'_'//trim(stp)//'.okc'
	    outfil=adjustl(outfil)
	    open(uijz(j),file=trim(outfil),status='unknown',err=10)
	    if(deepz>0) then
	      outrow=imax*(nzs+3)
	    else
	      outrow=imax*(nzs+2)
	    end if
	    write(uijz(j),*,err=11) nvar, outrow, 12 ! last entry on this line not used in all XMDV implemetations, but needed for VisIt to read database correctly RLB 2/21/2012
	    write(uijz(j),*,err=11) 'x'
	    write(uijz(j),*,err=11) 'y'
	    write(uijz(j),*,err=11) 'z'
	    write(uijz(j),*,err=11) 'p'
	    write(uijz(j),*,err=11) 'th'
	    write(uijz(j),fmt='(2(g18.9,1x),i3)',err=11) xmin,xmax,10
	    write(uijz(j),fmt='(2(g18.9,1x),i3)',err=11) ymin,ymax,10
	    write(uijz(j),fmt='(2(g15.5,1x),i3)',err=11) zlo,zhi,10
	    write(uijz(j),fmt='(2(g11.4,1x),i3)',err=11) plo,pmax,10
	    write(uijz(j),fmt='(2(g11.4,1x),i3)',err=11) th_min,th_max,10
	  end do
	end if
!
	return
  10	continue
  	write(*,*) 'Error opening output file ',outfil
  	write(ulog,*) 'Error opening output file ',outfil
  	close (ulog)
  	stop '10 in subroutine prpijz()'
  11	continue
  	write(*,*) 'Error writing output file ',outfil
  	write(ulog,*) 'Error writing output file ',outfil
  	close (ulog)
  	stop '11 in subroutine prpijz()'
	end subroutine prpijz
