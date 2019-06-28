c !  
c !  Reads a list of integer values that represent a 2-D matrix into 
c !  two 1-D arrays.  One array is a pointer that tracks the starting
c !  location of each row in the matrix.  This scheme is advantageous 
c !  for storing sparse "sawtooth" arrays. Such arrays have rows of 
c !  variable length and all of the non-zero values are at the left 
c !  ends of the rows. 
c !  by Rex L. Baum, spring 2001 
c !  
c ! Implementation:
c !  	open (u,file=infil,status='old')
c !  	call irdswm(len,jmax,u,test,x,u2,ctr)
c !  	close(u)
c !  
c ! 
  	subroutine irdswm(len,jmax,u,test,x,ctr,ulog)
	implicit none
	integer j,k,jmax,u,val,test,len,x(len),ulog
	integer ctr(jmax+1)
	j=0
	ctr(1)=0
	k=0
  100	continue
	read (u,*,end=110) val
c !  Check for marker value at end of a row
	if (val.eq.test) then
	read (u,*,end=110) j
	ctr(j)=k+1
	if (j.gt.jmax) go to 200
	go to 100
	end if
	k=k+1
	x(k)=val
	go to 100
  110	continue
 	ctr(jmax+1)=k+1
  	return
  200	continue
 	write(*,*) 'Array bound <jmax> exceeded in subroutine irdswm'
 	write(ulog,*) 'Array bound <jmax> exceeded in subroutine irdswm'
  	write(*,*) 'Press RETURN to exit'
  	read*
	stop '-200 in irdswm()'
	end