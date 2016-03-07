c !  
c !  subroutine to identify cells downslope from each grid cell
c !   and compute weighting factors for partitioning runoff among
c !   the cells
c !  
c !  Rex L. Baum, USGS, 4 March 2002, Latest revision 13 Mar 2013, RLB 
c !  
	subroutine slofac(nrow,ncol,rc,z,celsiz,nodat,nodata,cel,
     +	u,dscfil,pwr,next,ordr,u2,wffil,dsctr,dir,ridge,spars)
	implicit none
	integer k0,k,m0,m,i,j,l,u,rc,u2,nodata
	integer nrow,ncol,cel(ncol,nrow),count
	integer num(9),next(rc),ctr,ordr(rc)
	integer lh(9),rh(9),ll,lr,lmax,dsctr,dir(ncol,nrow)
	integer ridge(rc),spars
	real z(ncol,nrow),pwr,steep
	real pi,pio4,dzdx1,dxdx1,dydx1,dzdx2,dxdx2,dydx2
	real r,bs,as,dipdir,nodats,test
	double precision celsiz,nodat,diag,slope(9)
	double precision sum1,sum2,w(9),a,b
	character dscfil*255, wffil*255, outfil*255
	character scratch*6
	data lh/4,1,2,7,5,3,8,9,6/
	data rh/2,3,6,1,5,9,4,7,8/
	diag=celsiz*dsqrt(2.D0)
	pi=3.1415926535
	pio4=45.*pi/180.
	nodats=nodat
c ! in general, dxdx1=cos(az_x1); dydx1=sin(az_x1); dxdx2=cos(az_x2);dydx2=sin(az_x2)
c ! if x1 is the direction of maximum slope (d8) for square grids, they simplify to 
	dxdx1=1.
	dydx1=0.
	dxdx2=sqrt(2.)/2.
	dydx2=sqrt(2.)/2.
	r=dxdx2/dxdx1
c ! open output files
	outfil=dscfil
	outfil=adjustl(outfil)
  	open (u,file=trim(outfil),err=300)
	outfil=wffil
	outfil=adjustl(outfil)
  	open (u2,file=trim(outfil),err=300)
c !  search all cells for slope direction
	dsctr=0
	do 100, i=1,nrow
	  do 90, j=1,ncol
	  test=abs(z(j,i)-nodats)
	  if (test.ge.0.1) then
	  write (u,*) nodata
	  write (u,*) cel(j,i)
	  write (u2,*) nodata
	  write (u2,*) cel(j,i)
	  count=0
	  do 50, k0=1,3
	  do 40, m0=1,3
c !  start at upper left corner and go around the cell left to right and
c ! top to bottom
	  k=k0-2
	  m=m0-2
	  count=count+1
c !  set slopes equal to 0 before computing to aid in sorting later
	  slope(count)=0.D0
	  num(count)=0
	  if (k+i.ge.1 .and. k+i.le.nrow)then
	   if (m+j.ge.1 .and. m+j.le.ncol)then
c !  if adjacent cell is a nodata cell, set slope to a positive value 
	      test=abs(z(j+m,i+k)-nodats)
	      if(test.lt.0.1) then
	        slope(count)=10.D0
	         go to 40
	      end if
	      num(count)=cel(j+m,i+k)
	    if(k.eq.0 .or.m.eq.0) then
c !  compute slopes so that downslope is negative
	     slope(count)=(z(j+m,i+k)-z(j,i))/celsiz
	     else
	     slope(count)=(z(j+m,i+k)-z(j,i))/diag
	    end if
	   end if
	  end if
   40	continue
   50	continue
c !  compute weighting factors 
c ! case 0, pwr is out of range, so use default, which is 
c ! routing the flow down the steepest cells (D8 method)
	if (pwr.gt.20.0) then
c ! next line added to correct previous omission, 12 Nov 2010 RLB	
	  write (u,*) next(cel(j,i))
	  write (u2,1020) 1.0
	  dsctr=dsctr+1
	  go to 200
	end if

c ! case 1--uniform distribution among downslope cells
	if (pwr.eq.0.) then 
	  sum1=0.
   	  do 145, l=1,9
	  w(l)=0.
	  if(slope(l).lt.0) w(l)=1.
	  sum1=sum1+w(l)
  145	  continue
	  go to 170
	end if
c ! case 2-- distribution among downslope cells proportional to slope
	if (pwr.eq.1.) then 
	  sum1=0.
   	  do 155, l=1,9
	  w(l)=0.
	  if(slope(l).lt.0) w(l)=slope(l)
	  sum1=sum1+w(l)
  155	  continue
	  go to 170
	end if
c ! case 3--distribution according to power-law for slope
c ! increasing power toward infinity tends to diverting all
c ! water down the steepest path

	if (pwr.gt.0.0 .and. pwr.le.20.0) then
	  sum1=0.
   	  do 165, l=1,9
	  w(l)=0.D0
	  if(slope(l).lt.0) then
	    a=abs(slope(l))
	    w(l)=a**pwr
	  end if
	  sum1=sum1+w(l)
  165	  continue
	  go to 170
  	end if
c ! case 4--distribution of flow according to a variation 
c ! on Tarboton's D-infinity method, 
c ! as written works only for square grids
  	if (pwr.lt.0.0) then
	  sum1=0.
	  steep=0.
   	  do 166, l=1,9
	    w(l)=0.D0
  166	  continue
	lmax=dir(j,i)
	steep=slope(lmax)
c ! flat areas 
  	  if(steep.ge.0) then
	    write (u,*) next(cel(j,i))
	    write (u2,1020) 1.0
	    dsctr=dsctr+1
	    go to 200
	  end if
c ! sloping areas	  
	  lr=rh(lmax)
	  ll=lh(lmax)
	    a=slope(ll)
	    b=slope(lr)
	  if (a.eq.b) then
	     write (u,*) next(cel(j,i))
	     write (u2,1020) 1.0
	     dsctr=dsctr+1
	     go to 200
	  end if
	  if (a.lt.0. .and. b.lt.0.)then
	    if(a.gt.b) then
	      dzdx1=slope(lmax)
	      dzdx2=b
	      bs=(dzdx2-r*dzdx1)/dydx2
	      as=(dzdx2-bs*dydx2)/dxdx2
	      dipdir=atan(bs/as)
c ! route indeterminite cases down steepest path	      
	      if (dipdir.lt.0. .or. dipdir.gt.pio4) then
	        write (u,*) next(cel(j,i))
	        write (u2,1020) 1.0
	        dsctr=dsctr+1
	        go to 200
	      else
	        w(lmax)=pio4-dipdir
	        w(lr)=dipdir
	        sum1=pio4
		 if (w(lr).gt.w(lmax)) then
	           w(lmax)=dipdir
	           w(lr)=pio4-dipdir
		 end if
	        go to 170
	      end if
	     else 
	      dzdx1=slope(lmax)
	      dzdx2=a
	      bs=(dzdx2-r*dzdx1)/dydx2
	      as=(dzdx2-bs*dydx2)/dxdx2
	      dipdir=atan(bs/as)
	      if (dipdir.lt.0. .or. dipdir.gt.pio4) then
	        write (u,*) next(cel(j,i))
	        write (u2,1020) 1.0
	        dsctr=dsctr+1
	        go to 200
	      else
	        w(lmax)=pio4-dipdir
	        w(ll)=dipdir
		 if (w(ll).gt.w(lmax)) then
	           w(lmax)=dipdir
	           w(ll)=pio4-dipdir
		 end if
	        sum1=pio4
	        go to 170
	      end if
	    end if
	  end if
c ! a>0 and b>0	  
	      write (u,*) next(cel(j,i))
	      write (u2,1020) 1.0
	      dsctr=dsctr+1
	      go to 200
	end if 
  170	continue  
c ! attempt to find ridge crests ....................  
          if(sum1 .ge. spars) ridge(cel(j,i))=1
c ! round weighting factors to three decimal places if (a.lt.b) then
c ! and fix total to equal 1.  
	ctr=0
	sum2=0.
	do 175, l=1,9
	  w(l)=w(l)/sum1
	  if (w(l).lt.0.001) w(l)=0.
	  if (w(l).gt.0.) then
	    ctr=ctr+1
	    write (scratch,1020) w(l)
	    read (scratch,*) w(l)
	  end if
	  sum2=sum2+w(l)
  175	continue
  	do 180 l=1,9
	if (num(l).eq.next(cel(j,i))) then
	  w(l)=w(l)+(1-sum2)
	end if
  180	continue
c ! if steepest slope is zero, diverts flow to the next
c ! slope in the steepest downslope path
  	if (ctr.eq.0) then
	    write (u,*) next(cel(j,i))
	    write (u2,1020) 1.0
	    dsctr=dsctr+1
	    go to 200
	end if	
	do 185, l=1,9
	  if (w(l).gt.0) then
	    write (u,*) num(l)
	    write (u2,1020) w(l)
	    dsctr=dsctr+1
	  end if
  185	continue
  200	continue
	end if
   90	continue
  100	continue
     	close(u)
	return
c ! error reporting	
  300	continue
   	write (*,*) '*** Error opening output file ***'
	write (*,*) '--> ',outfil
	write (*,*) 'Check file path and status'
                write(*,*) 'Press RETURN to exit'
                read*
	stop
 1015	format(f6.0)
 1020	format(f6.3)
	end
	
