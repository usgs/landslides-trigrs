!  Implementation of Iverson's (2000) method of computing pore pressure
!  for rain infiltration.
!  by W.Z. Savage, spring 2001, with modifications by R.L. Baum
!  (both) USGS, Latest revsion 29 Jan 2013, RLB 
!  
	subroutine ivestp(u1,rikf,ulog,i,rf)
	use grids; use input_vars
	use model_vars
	implicit none
	integer:: i,j,jf,u1,ulog,nmx ! ,nccs
	integer:: n,nn
	real:: rikf(nts+1),finf 
	real (double):: derfc,a1,b1,ff,zns,zinc,z,t0,znew
	real (double):: zstar,tstar,x1,x2,x3,x4
	real (double):: rf1,rf2,rf3,rf4,rfa,rfb,rf(nzs+1)
	real (double):: fs,rslo,rphi,fmn,ptest,pmn,dhat ! pmn added 4/15/2010, dhat added 8 Jan 2013
	real (double):: captstar1,captstar2,tdif1,tdif2
	real (double):: dusz1,ddg2rad,newdep
!	logical:: lcvs 
	pi=3.141592653589793
	ddg2rad=pi/180.D0
!  maximum value of Factor of Safety	
 	finf=10.
 	nmx=0
! 	nmn=1+mmax ! initialization moved to main program 07 Jan 2013 RLB 
 	rslo=slo(i) 
	rphi=phi(zo(i))
	a1=sin(rslo)
	b1=cos(rslo)
	dhat=4.*dif(zo(i))/(b1*b1) ! added 8 Jan 2013, RLB
	p0zmx=0. ! Added 6 May 2013, RLB 
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
	dusz1=0.
	ptran=0. ! added 17Aug2009 RLB
	temporal_loop: do n=1,nts+1
	  t0=tcap(n)
	  fmn=1.e25
	  jf=jsav(n)
!	  if (flag<=-1 .and. flag>=-3) then !12/14/10 added flag>=-3
!	    write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
!	    'cell ',i,rslo/ddg2rad,n,t0 ! added label 4/22/2010 RLB
!	  end if
	  z=zmin
	  Z_loop: do j=1,nzs+1
	    znew=z
	    if(znew < 1.0e-30) znew =1.0e-30
	    if (abs(a1)>1.e-5) then
	      fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
	    else
	      fc(j)=0.d0
	    end if
	    pzero(j)=beta*(z-depth(i))
!	    if(z<dusz1 .or. rikf(n)==0.) then ! Revised 23 Oct 2013, RLB
	    if(z<dusz1) then
	      rf(j)=0.0
	    else
	      rf(j)=0.0
	      if (abs(z)>0.) then ! Formulas in next 2 lines apply only for z>0, test added 21 Feb 2013, RLB
	        zstar=z**2/dhat ! formula simplified 8 Jan 2013, RLB
	        tstar=t0/zstar
	      end if
	      temporal_loop_1: do nn=1,nper
	        if(z==0.) then ! exact formula added for case of z=0, 8 Jan 2013, RLB
	          tdif1=t0-capt(nn)
	          if(tdif1>0.) then
	            rfa=sqrt(tdif1*dhat/pi)
	          else
	            rfa=0.0
	          end if
	          tdif2=t0-capt(nn+1)
	          if(tdif2>0.) then
	            rfb=sqrt(tdif2*dhat/pi)
	          else
	            rfb=0.0
	          end if
	        else ! z>0
	          captstar1=capt(nn)/zstar 
	          tdif1=tstar-captstar1
	          if(tdif1 > 0.0) then 
	            x1=1./tdif1
	            x2=1./(sqrt(tdif1))
 	            rf1=sqrt(1./(x1*pi))*exp(-x1)
  	            rf2=derfc(x2)
	            rfa=rf1-rf2
	          else
	            rfa=0.0
	          end if
	          captstar2=capt(nn+1)/zstar
	          tdif2=tstar-captstar2
	          if(tdif2 > 0.0) then	
	            x3=1./tdif2
	            x4=1./(sqrt(tdif2))
	            rf3=sqrt(1./(x3*pi))*exp(-x3)
	            rf4=derfc(x4)
	            rfb=rf3-rf4
	          else
	            rfb=0.0
 	          endif
 	        endif
 	        rf(j)=rf(j)+rik(i+(nn-1)*imax)*(rfa-rfb)
   	        if(rfa==0.0 .and. rfb==0.0) exit ! skip unnecessary cycles, RLB, 2/19/2015
! 	        if (z==0) write(*,*) 'i, j, rfa, rfb, rf(j), tdif1, tdif2',  i, j, rfa, rfb, rf(j), tdif1, tdif2
  	      end do temporal_loop_1
	    end if
	    bline(j)=z*beta
	    if(abs(rf(j))>0.0) then ! added 17AUG2009 RLB 
              ptran(j)=z*rf(j) ! correction, added "z*" 4/27/2012, RLB
              if(z==0) ptran(j)=rf(j) ! formula for ptran at z=0 not normalized, 9 Jan 2013, RLB 
	    end if
 	    p(j)=pzero(j)+ptran(j)
	    ptest=p(j)-bline(j)
	    if(ptest > 0.0) then
	      p(j)=bline(j)
	    end if
	    z=z+zinc
   	  end do   Z_loop
   	  if(n==1) p0zmx=p(nzs+1) ! Added 6 May 2013, RLB 
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
! Compute factor of safety & save results  	  
   	  z=zmin
	  Z_FS_loop: do j=1,nzs+1 
	    if (abs(a1)>1.e-5 .and. z>1.e-30) then
	      if(lpge0 .and. p(j)<0.) then !option added 4/15/2010
	        fw(j)=0.d0
	      else
	        fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
	      end if
	    else
	      fw(j)=0.d0
	    end if
	    fs=ff+fw(j)+fc(j)
! frictional strength cannot be less than zero 	
	    if ((ff+fw(j))<0.) fs=fc(j)
 	    if (fs>finf) fs=finf
 	    if (z<=1.e-02) fs=finf 
! 	    if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
! 	    if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j)&
!            & ,pzero(j),ptran(j),bline(j),fs
! 	    if (flag==-3 .and. unsat0) then ! revised 12/23/2010
! 	      write(u1,'(6(g12.5,1x):)') z,p(j),fs,1. ! added 4/14/2010 RLB
! 	    else if (flag==-3) then
! 	      write(u1,'(6(g12.5,1x):)') z,p(j),fs
! 	    end if
	    if(jf>0) then
	      if (fs<fmn) then
	        fmn=fs
	        zfmin(i+(jf-1)*imax)=z
	        pmn=p(j) ! revised 4/15/2010
	      end if
! Store pressure head and related output in 3-d arrays. Added 17Nov2014, RLB
 	        if(flag<0 .or. outp(1)) then
                  p3d(i+(jf-1)*imax,j)=p(j)
                  newdep3d(i+(jf-1)*imax)=newdep
                  dh3d(i+(jf-1)*imax)=0.d0
                end if 
 	      if(flag==-1) fs3d(i+(jf-1)*imax,j)=fs
 	      if(flag==-2) then
 	        fs3d(i+(jf-1)*imax,j)=fs
 	        ptran3d(i+(jf-1)*imax,j)=ptran(j)
 	        pzero3d(i,j)=pzero(j)
              end if
 	      if(flag==-3) then
 	        fs3d(i+(jf-1)*imax,j)=fs
 	        th3d(i+(jf-1)*imax,j)=ths(zo(i))
              end if
              if(flag<=-4 .or. outp(1)) th3d(i+(jf-1)*imax,j)=ths(zo(i))
	    end if
	    z=z+zinc
   	  end do Z_FS_loop
	  if (jf>0) then  ! revised 4/29/2010 to include pmin() RLB
	    fsmin(i+(jf-1)*imax)=fmn
	    if(fmn==finf) then ! Added 30 Jan 2013, RLB 
	      pmn=p(nzs+1)
	      zfmin(i+(jf-1)*imax)=zmax(i)
            end if 
	    if(lpge0 .and. pmn<0.) then !option added 4/15/2010
	      pmin(i+(jf-1)*imax)=0.
	    else
	      pmin(i+(jf-1)*imax)=pmn
	    end if
!	    if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
!  	      dcf=0. ! this line used only for saturated infiltration model
!  	      chi=1.d0 !added 12/23/2010
!	      if (flag>=-6) call svijz(i,jf,0.d0,newdep,ulog)
!	      if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,0.d0,newdep,ulog)  ! Added 2/10/2012
!	    end if
	  end if
  	end do temporal_loop
	return
	end subroutine ivestp
