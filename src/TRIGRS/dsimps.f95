	subroutine dsimps(n,h,f,intf)
! by Rex L. Baum 7 July 2004,	Latest revision 29 Jan 2013, RLB
! uses Simpson's 3-point rule to numerically integrate a function at n successive points
! on an interval that has been subdivided into 2n equal segments.
	integer, parameter:: double=kind(1d0)
	integer:: i,n
	real (double):: f(2*n+1)
	real (double):: h,intf(n+1),sumf,fl,fmid,fr ! n changed to n+1, 29 Jan 2013, RLB
	sumf=0.d0
	do i=1,n
	  fl=f(2*i-1)
	  fmid=f(2*i)
	  fr=f(2*i+1)
	  sumf=sumf+fl+4.d0*fmid+fr
	  intf(i+1)=sumf*h/3.d0 ! value of integral over the interval f(1) and f(2*i+1)
	end do
	return
	end subroutine dsimps
	