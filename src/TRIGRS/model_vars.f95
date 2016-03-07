! By Rex L. Baum, 1 April 2004
module model_vars
	integer,parameter:: double=kind(1d0)
	integer:: nts,kper,nmax1,nmax2,nmn,nmin
	integer,allocatable:: jsav(:), ix(:),jy(:) ! added ix() & jy() 4/21/2010 
	real:: test1,dg2rad
	real,allocatable:: q(:),qb(:)
	real (double):: eps,tmin,tmax,ts,qt,tns,beta,qmax,tinc
	real (double):: test,nodat,sumex,dusz,dcf,vf0,p0zmx ! Added p0zmx 6 May 2013, RLB
	real (double):: celsiz,param(6),parami(6),ti,tis,pi,smt,lard,xllc,yllc,zmn(1),zmx(1) ! added xllc & yllc 12/24/2010
	real (double),allocatable:: p(:),ptran(:),pzero(:),bline(:),chi(:) ! added chi 12/22/2010 RLB
	real (double),allocatable:: r(:),fc(:),fw(:),thz(:),kz(:),tcap(:),tinc_sat(:) ! tinc_sat added 25 June 2013, RLB 
	real (double),allocatable:: trz(:),uwsp(:),gs(:),qtime(:),qts(:) ! Added qts 21 Feb 2013, RLB 
	real (double),allocatable:: p3d(:,:),pzero3d(:,:),ptran3d(:,:),fs3d(:,:),th3d(:,:),dh3d(:),newdep3d(:) ! Added 17Nov2014, RLB
end module model_vars

