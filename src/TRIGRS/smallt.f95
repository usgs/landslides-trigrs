	subroutine smallt(tstar,alfa,qa,qb,d,&
        &qlb,ths,thr,qt,psih,ck,th,z)
! By Rex L. Baum and W.Z. Savage, USGS, October 1, 2004        
	integer,parameter:: double=kind(1d0)
     	real (double), intent(in):: tstar,alfa,qa,d,qlb,z ! added intent 12/22/2010, RLB ,ths,thr,qb
     	real, intent(in):: ths,thr,qb ! corrected to single 1/21/2011, 2/13/2013 RLB
     	real (double), intent(out):: qt,psih,ck,th ! added intent 12/22/2010, RLB
     	real (double):: ck0,psi0,th0,f0,ft,gt,arg1,arg2 
     	real (double):: f1a,f1b,f1c,f1,f2a,f2b,f2c,f2
     	real (double):: gt1,gt2,derfc
        ck0=qa-(qa-qlb)*exp(-alfa*(d-z))
 	psi0=(log(ck0))/alfa
     	th0=thr+(ths-thr)*exp(alfa*psi0)
   	f0=2.*(qb-qa)*exp(alfa*z/2.)
    	ft=0.0
	gt=0.0
	arg1=(alfa*z)/(2.*sqrt(tstar))
	f1a=derfc(arg1)
	f1b=exp(arg1*sqrt(tstar)+tstar/4.)
	f1c=derfc(arg1+sqrt(tstar/4.))
	f1=f1a-f1b*f1c
	arg2=alfa*(2.*d-z)/(2.*sqrt(tstar))
	f2a=derfc(arg2)
	f2b=exp(arg2*sqrt(tstar)+tstar/4.)
	f2c=derfc(arg2+sqrt(tstar/4.))
	f2=-f2a+f2b*f2c
	ft=f1+f2 
	gt1=f1b*f1c
	gt2=f2b*f2c
	gt=gt1-gt2
        ck=ck0+f0*ft
	if(ck < 0.0) ck=-ck
   	psih=log(ck)/alfa
	th=thr+(ths-thr)*exp(alfa*psih)
	qt=f0*gt/2.0
	return
	end subroutine smallt
