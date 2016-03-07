	subroutine dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
c ! By Rex L. Baum, May 14, 2004
c ! Uses bisection method as described in Press and others, 1986, p. 246-247.
	implicit none
	integer imax,m,mmax
	double precision a,b,x0,eps
	double precision rtb,fm,dx
	imax=0
  	do 10,m=1,mmax
  	  dx=dx/2.
  	  x0=dx+rtb
  	  fm=sin(a*x0)+b*x0*cos(a*x0)
  	  if(fm.le.0.) rtb=x0
  	  if(abs(dx).le.eps .or. fm.eq.0.) then
  	    imax=m
  	    return
  	  end if 
   10  	continue 
   	return
   	end
