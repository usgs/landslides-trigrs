	subroutine sindex(n,ra,ndx)
c ! uses heapsort algorithm	
c ! based on section 8.3 of Numerical Recipes, 
c ! Indexes a real array "ra" of length n,
c ! i.e. outputs the array ndx such that
c ! ra(ndx(j)) is in ascending order for j=1,2, ...,n.  
c ! The input quantities n and ra are not changed.
	implicit none
	integer n, ndx(n), j,l,ir,ndxt,i
	real ra(n),q
	
	do 11, j=1,n
	ndx(j)=j
   11	continue
	l=n/2+1
	ir=n
   10	continue
   	if (l.gt.1) then
	  l=l-1
	  ndxt=ndx(l)
	  q=ra(ndxt)
	else
	  ndxt=ndx(ir)
	  q=ra(ndxt)
	  ndx(ir)=ndx(1)
	  ir=ir-1
	  if(ir.eq.1) then
	  	ndx(1)=ndxt
		return
	  endif
	endif
	i=l
	j=2*l
   20	if(j.le.ir)then
   	  if(j.lt.ir) then
	     if(ra(ndx(j)).lt.ra(ndx(j+1)))j=j+1
	   end if
	   if(q.lt.ra(ndx(j)))then
	     ndx(i)=ndx(j)
	     i=j
	     j=j+j
	   else
	     ndx(i)=ndxt
	     go to 10
	   endif
	go to 20
	end if
	ndx(i)=ndxt
	go to 10
	end
