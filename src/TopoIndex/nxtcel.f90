!  
!  subroutine to identify subjacent cells for a DEM grid
!  uses a grid of slope directions for input
!  
!  Rex L. Baum, USGS, 15 February 2002, 
!  29 Sept. 2005 converted to Fortran 90
! 
	subroutine nxtcel(nrow,ncol,rc,prm,u1,nodat,nodata,&
        &z,dir,cell,list,cels,u2)
	implicit none
	integer nrow,ncol,prm,u1,i,j,kount,inew,jnew,a,b,k,rc
	integer m(9),n(9),cell(ncol,nrow),dir(ncol,nrow),list(prm)
	integer cels(rc),u2,nodata,mmcnt
	real z(ncol,nrow),nodats,test
	double precision nodat
	data m/-1,-1,-1,0,0,0,1,1,1/
	data n/-1,0,1,-1,0,1,-1,0,1/
	kount=0
	mmcnt=0
	nodats=nodat
	write(*,*) 'nxtcel() nodata (integer,floating)= ',nodata,nodats
	do i=1,nrow
	  do j=1,ncol
	   test=abs(z(j,i)-nodats)
	   if (test >= 0.1) then
	     kount=kount+1
	     cell(j,i)=kount
	   end if
   	  end do
  	end do
	write(*,*) 'Identifying downslope cells and grid mismatches'
	write(u2,*) ''
	write(u2,*) '         Listing of grid mismatches'
	write(u2,*) 'Mismatch counter, Row, Column, Direction code'
	Row_traverse: do i=1,nrow
	  Column_traverse: do j=1,ncol
	   test=abs(z(j,i)-nodats)
	   if (test < 0.1) then
	     cycle Column_traverse 
   	   end if
   	   kount=0
	   inew=i
	   jnew=j
  	   Nextcell: do
  	     if(dir(jnew,inew)==nodata) then
	       mmcnt=mmcnt+1
	       write (u2,*) mmcnt, jnew, inew, dir(jnew,inew)
	       cycle Column_traverse 
	     end if
	     if (dir(jnew,inew).ne.5.or.dir(jnew,inew).ne.0) then
	       if(inew>=1.and.inew<=nrow.and.&
     	       &jnew>=1.and.jnew<=ncol) then
	         kount=kount+1
!  compute indices of next cell
	         a=inew+m(dir(jnew,inew))
	         b=jnew+n(dir(jnew,inew))
	         if (a==inew.and.b==jnew) then
	           exit nextcell 
	         end if
	         if(a>=1.and.a<=nrow.and.b>=1.and.b<=ncol) then
	           list(kount)=cell(b,a)
	           if(list(kount)==0) list(kount)=cell(jnew,inew)
	           test=abs(z(jnew,inew)-nodat)
	           if (test<=0.1) exit nextcell 
	           if (kount==2) exit nextcell 
 	           cycle nextcell 
	         else
	           exit nextcell 
	         end if
	       end if
	     end if
	   end do nextcell
	      if(u1>0) write(u1,*) nodata
	      if(u1>0) write(u1,*) cell(j,i)
	      if (kount>1) then
	        cels(cell(j,i))=list(1)
	        do k=1,kount-1
	         if(u1>0) write(u1,*) list(k)
  	        end do
   	      else
	       if(u1>0) write(u1,*) cell(j,i)
	       cels(cell(j,i))=cell(j,i)
	      end if
	  end do Column_traverse 
  	end do Row_traverse
  	if (mmcnt>0) then
	  write(*,*) mmcnt, ' mismatched cells'
	  write(*,*) 'Reconcile elevation grid and direction grid'
	  write(*,*) 'before attempting to use TopoIndex'
	  write(u2,*) mmcnt, ' mismatched cells'
	  write(u2,*) 'Reconcile elevation grid and direction grid'
	  write(u2,*) 'before attempting to use TopoIndex'
                  write(*,*) 'Press RETURN to exit'
                  read*
	  close (u2)
	  stop
	end if
 	write (u2,*) mmcnt,',  --,   --,  --'
 	write (*,*) 'No grid mismatch found!'
 	write (u2,*) 'No grid mismatch found!'
 	write (u2,*) 'Subroutine nxtcel completed normally'
	return
	end
