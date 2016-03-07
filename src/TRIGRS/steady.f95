 	subroutine steady(sumex,ulog,imx1)
 	use grids; use input_vars
 	implicit none
 	integer:: i,acnt,ulog,imx1
	integer,parameter:: double=kind(1d0)
	real (double):: sumex,b,rslo
!  By Rex L. Baum, 1 April 2004
!  Compute initial estimate of Isteady/Ks, rikzero(), and test values.
!  Isteady must be < Ks.  If cos(slo) *cos(slo)<(Isteady/Ks), then an 
!  inverted water table results.
	write(*,*) 'Testing and adjusting steady infiltration rates'
	acnt=0
	do i=1,imx1
	  if(ks(zo(i))==0.) then ! prevent division by zero errors
	    rikzero(i)=1.
	  else
	    rikzero(i)=rizero(i)/ks(zo(i))
	  end if
	end do  
	sumex=0
	do i=1,imx1
	  rslo=slo(i)
	  b=cos(rslo)*cos(rslo)
	  if(rikzero(i)>=b) then
	    rikzero(i)=cos(rslo) 
	    acnt=acnt+1
	    write(ulog,*) 'Adjusted steady infiltration rate, cell ',i
	    write(ulog,*) 'Corrected rizero/Ks = ',rikzero(i)
	    write(ulog,*) 'Original rizero, ks',rizero(i),ks(zo(i))
	    write(ulog,*) 'Set pore presssures to zero'
	  end if
	  if (depth(i)==0 .and. rizero(i)<0) then
	  sumex=sumex-rizero(i)
	  end if
	end do
	write(*,*) 'Adjusted steady infiltration rate at '&
        &,acnt,' cells'
	write(ulog,*) 'Adjusted steady infiltration rate at '&
        &,acnt,' cells'
	return
	end subroutine steady