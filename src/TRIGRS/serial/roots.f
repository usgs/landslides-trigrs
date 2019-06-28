	subroutine roots(nmax,r,d,beta,eps,pi)
c ! by Rex L. Baum and W.Z. Savage, May 2004
c! Latest revision 28 Jan 2013, RLB 	
	implicit none
	integer i,n,nmax,imax,mmax
c !	real beta promoted to double precision 12/6/2006, df1 & dfsgn added 29 Jan 2013, RLB 
	double precision a,b,x0,f,df,x1,d,df1,dfsgn
	double precision r(*),pi,eps,beta
	double precision rlb,rub,rtb,fl,fu,dx
	a=beta*d
	b=2.0
	x0=eps
	do 15 n=1,nmax
c !  set upper and lower bound for each interval
	  rlb=(pi/a)*float(n-1)+eps
	  rub=(pi/a)*float(n)
  	  fl= sin(a*rlb)+b*rlb*cos(a*rlb)
  	  fu= sin(a*rub)+b*rub*cos(a*rub)
  	  if(fl.lt.0) then
  	    rtb=rlb
  	    dx=rub-rlb
  	  else
  	    rtb=rub
  	    dx=rlb-rub
  	  end if
c ! find initial estimate of roots by bisection	
	  mmax=5
	  call dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
	  if(imax.ge.1) go to 11
  	  do 10 i=1,100
  	    f= sin(a*x0)+b*x0*cos(a*x0)
c !	    df=(a-b)*cos(a*x0)-a*b*x0*sin(a*x0) !formula appears to contain sign error 28 Jan 2013, RLB
	    df=(a+b)*cos(a*x0)-a*b*x0*sin(a*x0)
c ! Next lines added to prevent division by zero, 28-30 Jan 2013, RLB. 	   
            df1=abs(df) 
	    if(df1 .lt. eps) then
              dfsgn=1.d0
              if(df .lt. 0.d0) dfsgn=-1.d0
              df=eps*dfsgn
            end if
	    x1=x0-f/df
	    if(abs((x1-x0)/x0) .lt. eps) then
	      if(x1.gt.rlb .and. x1.lt.rub) then
	        imax=i
	        x0=x1
  	        goto 11
  	      else
c ! if Newton-Raphson method converges on the wrong root, then
c ! revert to bisection
	        mmax=200
	        call dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
  	        imax=imax+5
	        if(imax.ge.6) go to 11
	      end if
	    else
	      x0=x1
  	    end if
  10	  continue
c ! if Newton-Raphson method did not converge, then
c ! revert to bisection
  	    mmax=200
	    call dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
  	    imax=imax+5
  11	  continue
  	  r(n)=x0
  	  x0=x0+pi/a
  15	continue
        return
  	end


