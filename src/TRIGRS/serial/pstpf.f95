!  diffusion model for pressure head from applied load of rising water table
!By R.L. Baum and W.Z. Savage, USGS, May 2004, latest revision 21 Jan 2011
subroutine pstpf(u1,dhwt,dwt,&
  & ulog,i,iper,rf,nccs,lcvs,nmx,t0,jf)
use grids; use input_vars
use model_vars
implicit none
integer:: i,j,jf,u1,nmx,ulog,n1,nccs,n,m,iper,jmark,nstp,jtop
real:: dwt !,chi 12/22/2010 made chi a 1-d array & moved to model_vars 
real (double):: finf,t1,t2,term1,term2,rn,znew,t0
real (double):: ferfc1,ferfc2,ferfc3,ferfc4,dhwt(nts+1),dh
real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
real (double):: ar1,ar2,ar3,ar4,fs,rf(nzs+1)
real (double):: rfa,rfb,derfc,ff,rslo,rphi,fmn,ptest
real (double):: tol,delt1,delt2,t1old,t2old
real (double):: ddg2rad,dusz0,dusz1,dlz,d1,zm
real (double):: tstar,flt1,flt2,uwt1,uwsum
logical:: lcvs 
pi=3.141592653589793
ddg2rad=pi/180.D0
!  maximum value of Factor of Safety
finf=10.
nmx=0
nmn=1+mmax
tol=1.e-06
dh=dhwt(iper+1)
rslo=slo(i)
rphi=phi(zo(i))
a1=sin(rslo)
b1=cos(rslo)
newdep=-9999. ! Added 4 Jan 2016, RLB
select case (flowdir) ! set value of beta (Iverson's beta line)
  case ('slope')
    beta=b1*b1
  case ('hydro')
    beta=1.d0
  case default
    beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
  end select
if(abs(b1-rikzero(i))<1.e-6) beta=0.d0 
if (abs(rslo)>1.e-5) then
  ff=tan(rphi)/tan(rslo)
else !  set factor of safety to fixed value for flat slopes   
  ff=finf
end if
zns=float(nzs)
zinc=(zmax(i)-zmin)/zns  
nstp=int((1./alp(zo(i))/zinc))
z=zmin
fmn=1.e25
! compute depth of unsaturated zone, 
if(unsat(zo(i))) then
  if(lps0) then
    dusz0=depth(i)-1./alp(zo(i)) ! use when psi0=-1/alp(zo(i))
  else
    dusz0=depth(i) ! use if psi0=0
  end if
    dusz1=dwt-1./alp(zo(i))
else
  dusz0=0.  
  dusz1=0.  
end if
dlz=zmax(i)-dusz0
zm=zmax(i)
uwsum=0. ! for computing depth-averaged unit weight
uwsp=uws(zo(i))
rf=0.0 
jmark=1 ! Added in case depth(i)==0., 28 Jan 2013, RLB 
Z_loop: do j=1,nzs+1
  znew=z
  if(z<depth(i))then 
    jmark=j
  end if
  if(znew < 1.0e-30) znew =1.0e-30
  pzero(j)=beta*(z-depth(i))
  if(z <= depth(i))then 
    pzero(j)=0.
    if(lps0) then ! psi0=-1/alpha
      if(llus) then ! compute water table rise
        if(z > dusz0-dh .and. z < depth(i)-dh) then  
          rf(j)=(-1./alp(zo(i))+(z-dusz0+dh))
        else if(z >= depth(i)-dh) then  
          rf(j)=beta*(z-depth(i)+dh)
        else 
          rf(j)=ptran(j)
        end if
      else ! static water table
        if(z > dusz0 .and. z < depth(i)) then  
          rf(j)=(-1./alp(zo(i))+(z-dusz0))
        else if(z >= depth(i)) then  
          rf(j)=beta*(z-depth(i))
        else
          rf(j)=ptran(j) 
        end if
      end if
    else ! psi0=0
      if(z >= dusz0-dh .and. llus) then ! water table rise
        rf(j)=beta*(z+dh-dusz0)
      else ! static water table
        rf(j)=ptran(j)
      end if
    end if
  else if(z>=zm .and. dlz<0.001) then
    rf(j)=beta*(z+dh-dusz0)      
  else if(llus) then  ! compute pressure rise below wt only if llus=.true.
    temporal_loop: do m=1,iper
      tdif1=t0-tcap(m)
      if(tdif1 > 0.0) then
        d1=dif(zo(i))/(b1*b1)
        t1=sqrt(d1*tdif1)
        tstar=t1/(dlz*dlz)
        if (t1<1.0e-29) t1=1.0e-29
        term1=0.0
        if(tstar>5.) then ! solution for later time...
          series_a_lt: do n=1,mmax
            rn=float(n)
            ar1=(2.d0*rn-1.d0)*pi*(z/dlz-1.d0)/2.d0
            ar2=(2.d0*rn-1.d0)*(2.d0*rn-1.d0)*pi*pi*tstar/4.d0
            flt1=cos((rn-1.d0)*pi)*exp(-ar2)*cos(ar1)/(2.d0*rn-1.d0)
! test for convergence of series       
            t1old=term1
            term1=term1+flt1
            delt1=abs(term1-t1old)
            if (term1/=0.) then
              delt1=delt1/abs(term1)
            end if
            n1=n
            if(delt1<=tol) exit
            if(abs(term1)<=tis .and. n>=4) then
              term1=0.d0
              delt1=0.d0
              exit
            end if
          end do series_a_lt
          term1=1.d0-term1*4.d0/pi
        else ! solution for early time
          series_a_et: do n=1,mmax
            rn=float(n)
            ar1=((2.*rn-1.)*dlz-(zm-z))/(2.*t1)
            ar2=((2.*rn-1.)*dlz+(zm-z))/(2.*t1)
            ferfc1=cos((rn+1.0)*pi)*derfc(ar1)
            ferfc2=cos((rn+1.0)*pi)*derfc(ar2)
! test for convergence of series       
            t1old=term1
            term1=term1+ferfc1+ferfc2
            delt1=abs(term1-t1old)
            if (term1/=0.) then
              delt1=delt1/abs(term1)
            end if
            n1=n
            if(delt1<=tol) exit
          end do series_a_et
        end if
        if(lcvs .and. delt1>tol) then
          nccs=nccs+1
          nv(i)=1
          lcvs=.false.
          write(*,*) 'SZ noncvg a',i,m,n1,t0,term1,t1old
          write(*,*) 'delt1,tol',delt1,tol
        end if
        if(n1>nmx) nmx=n1
        if(n1<nmn) nmn=n1
        rfa=term1
      else
        rfa=0.0
      end if
      tdif2=t0-tcap(m+1)
      if(tdif2 > 0.0) then
        d1=dif(zo(i))/(b1*b1)
        t2=sqrt(d1*tdif2)
        tstar=t2/(dlz*dlz)
        if (t2<1.0e-29) t2=1.0e-29
        term2=0.0 
        if(tstar>5.) then ! solution for later time...
          series_b_lt: do n=1,mmax
            rn=float(n) ! Added 14Sep2009, RLB, previously omitted in error.
            ar1=(2.d0*rn-1.d0)*pi*(z/dlz-1.d0)/2.d0
            ar2=(2.d0*rn-1.d0)*(2.d0*rn-1.d0)*pi*pi*tstar/4.d0
            flt2=cos((rn-1.d0)*pi)*exp(-ar2)*cos(ar1)/(2.d0*rn-1.d0)
! test for convergence of series      
            t2old=term2
            term2=term2+flt2
            delt2=abs(term2-t2old)
            if (term2/=0.) then
              delt2=delt2/abs(term2)
            end if
            n1=n
            if(n>10) write (*,*) 'b_lt: n,term2,flt1,tstar',n,term2,flt2,tstar
            if(delt2<=tol) exit
            if(abs(term2)<=tis .and. n>=4) then
              term2=0.d0
              delt2=0.d0
              exit
            end if
          end do series_b_lt
          term2=1.d0-term2*4.d0/pi
        else ! solution for early time
          series_b_et: do n=1,mmax 
            rn=float(n)
            ar3=((2.*rn-1.)*dlz-(zm-z))/(2.*t2)
            ar4=((2.*rn-1.)*dlz+(zm-z))/(2.*t2)
            ferfc3=cos((rn+1.0)*pi)*derfc(ar3)
            ferfc4=cos((rn+1.0)*pi)*derfc(ar4)
! test for convergence of series       
            t2old=term2
            term2=term2+ferfc3+ferfc4
            delt2=abs(term2-t2old)
            if (term2/=0.) then
              delt2=delt2/abs(term2)
            end if
            n1=n
            if(delt2<=tol) exit
          end do series_b_et
        end if
        if(lcvs .and. delt2>tol) then
          nccs=nccs+1
          nv(i)=1
          lcvs=.false.
          write(*,*)  'SZ noncvg b',i,m,n1,t0,term2,t2old
          write(*,*) 'delt2,tol',delt2,tol
        end if
        if(n1>nmx) nmx=n1
        if(n1<nmn) nmn=n1
        rfb=term2
      else
        rfb=0.0
      end if
      rf(j)=rf(j)+dh*beta*(rfa-rfb)
      if(rfa==0.0 .and. rfb==0.0) exit ! skip unnecessary cycles, RLB, 2/19/2015
    end do temporal_loop 
  end if
  bline(j)=z*beta
  ptran(j)=rf(j)
  p(j)=pzero(j)+ptran(j)
  ptest=p(j)-bline(j)
  if(ptest > 0.0) then
  p(j)=bline(j)
end if
! estimate partially saturated unit weight      
if(p(j)<-1.d0/alp(zo(i))) then
  uwt1=(gs(zo(i))*(1-ths(zo(i)))+thz(j))*uww
else
  uwt1=uws(zo(i))
end if
uwsum=uwsum+uwt1
uwsp(j)=uwsum/float(j)
z=z+zinc
end do Z_loop  
! find new height of rising water table in zones of upward seepage   
if(rikzero(i)<0.0) then
  zinc=(zmax(i)-zmin)/zns
  z=zmin
  newdep=0.0
  do j=1,nzs+1
    if(p(j)<0.0) newdep=z
    z=z+zinc
  end do
! adjust presures
  z=zmin
  do j=1,nzs+1
    if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
    if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
    z=z+zinc    
  end do
end if
! Smooth data near initial water table
if(llus .and. zm>depth(i)) then
  jtop=jmark-nstp
  if(jtop<1) jtop=1
  do j=jmark,jtop,-1 ! corrected 5/12/2014, RLB
    if(ptran(j)>ptran(j+1)) then 
      ptran(j)=ptran(j+1)
      p(j)=ptran(j)+pzero(j)
    end if
  end do
end if  
! Compute factor of safety & save results    
z=zmin
Z_FS_loop: do j=1,nzs+1    
  chi(j)=1.0 
! Approximate suction stress computation, moved 4/22/2010 RLB
  if(p(j)<0.) then ! for suction 2/25/2009 (formerly for suctions exceeding the air-entry value, -1/alpha)
    chi(j)=(thz(j)-thr(zo(i)))/(ths(zo(i))-thr(zo(i)))
  else
    chi(j)=1.0; thz(j)=ths(zo(i)) ! Revised 8/22/13 RLB 
  end if
  if (abs(a1)>1.e-5 .and. z>0.) then
    fw(j)=-(chi(j)*p(j)*uww*tan(rphi))/(uwsp(j)*z*a1*b1) ! same as saturated formula with factor Chi
    fc(j)=c(zo(i))/(uwsp(j)*z*a1*b1) ! uses depth averaged unit weight
  else
    fw(j)=0.d0
    fc(j)=0.d0
  end if
  fs=ff+fw(j)+fc(j)
! frictional strength cannot be less than zero 
  if ((ff+fw(j))<0.) fs=fc(j)
  if (fs>finf) fs=finf
  if (z<=1.e-02) fs=finf 
! Store pressure head and related output in 3-d arrays. Added 17Nov2014, moved 02Dec2014, RLB
  if (jf>0) then
    if(flag<0 .or. outp(1)) then
      p3d(i+(jf-1)*imax,j)=p(j)
      newdep3d(i+(jf-1)*imax)=newdep
      dh3d(i+(jf-1)*imax)=dh
    end if 
    if(flag==-1) fs3d(i+(jf-1)*imax,j)=fs
    if(flag==-2) then
      fs3d(i+(jf-1)*imax,j)=fs
      ptran3d(i+(jf-1)*imax,j)=ptran(j)
      pzero3d(i,j)=pzero(j)
    end if
    if(flag==-3) then
      fs3d(i+(jf-1)*imax,j)=fs
      th3d(i+(jf-1)*imax,j)=thz(j)
    end if
    if(flag<=-4 .or. outp(1)) th3d(i+(jf-1)*imax,j)=thz(j)
  end if
  if (fs<fmn) then
    fmn=fs
    if (jf>0) then
      zfmin(i+(jf-1)*imax)=z
      pmin(i+(jf-1)*imax)=p(j)
    end if
  end if
  z=z+zinc
end do Z_FS_loop
if (jf>0) then ! Revised 12/22/2010 RLB
  fsmin(i+(jf-1)*imax)=fmn
  if(fmn==finf) then ! Added 30 Jan 2013, RLB 
    pmin(i+(jf-1)*imax)=p(nzs+1)
    zfmin(i+(jf-1)*imax)=zmax(i)
  end if 
end if
return
end subroutine pstpf
