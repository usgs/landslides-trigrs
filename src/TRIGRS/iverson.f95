!  Implementation of Iverson's (2000) method of computing pore pressure
!  for rain infiltration.
!  by W.Z. Savage, spring 2001, with modifications by R.L. Baum
	subroutine iverson(imx1,u1,profil,ulog)
	use grids; use input_vars
	use model_vars
 	implicit none
	integer:: j,i,u1,ulog,imx1,n 
	real:: finf 
	real (double) :: derfc,a1,b1,ff,zns,zinc,z,znew
	real (double) :: zstar,tstar,x1,x2,x3,x4
	real (double) :: rf1,rf2,rf3,rf4,rfa,rfb,rf
	real (double) :: fs,rslo,rphi,fmn,ptest,pmn,dhat ! pmn added 4/15/2010, dhat added 8 Jan 2013
	real (double) :: newdep,captstar1,captstar2,tdif1,tdif2
	character (len=255) profil
	write(ulog,*) 'Starting saturated-zone'
	write(ulog,*) 'computations for infinite-depth'
	write(*,*) 'Starting saturated-zone'
	write(*,*) 'computations for infinite-depth'
	pi=3.141592653589793
	dg2rad=pi/180.D0
!  maximum value of Factor of Safety	
	finf=10.
!  loop steps through all grid cells
	write(*,*) 'Cells completed: '
	grid_loop: do i=1,imx1
	  rslo=slo(i) 
	  if(rslo<slomin .or. rslo>slomax .or. zmax(i)<=0.0001) then
	    fsmin(i)=finf+1.
	    zfmin(i)=zmax(i)
	    pmin(i)=0.
	    if (mod(i,2000)==0) write (*,fmt='(2x,i10,2x,a1)',advance="no") i,char(13) ! cells completed
	    cycle grid_loop
	  end if
!	  if (flag<=-1 .and. flag>=-3) write (u1,'(a5,i12,f6.1,i2,2x,g14.8)')&
!	  & 'cell ',i,rslo/dg2rad,1,t ! added label, t 4/22/2010, 12/14/10 added flag>=-3
	  rphi=phi(zo(i))
	  a1=sin(rslo)
	  b1=cos(rslo)
	  dhat=4.*dif(zo(i))/(b1*b1) ! added 8 Jan 2013, RLB
	  newdep=-9999. ! Added 4 Jan 2016, RLB
	  select case (flowdir) ! set value of beta (Iverson's beta line)
	    case ('slope')
	    beta=b1*b1
	    case ('hydro')
	    beta=1.d0
	    case default
	    beta=b1*b1-rikzero(i)
	  end select
	  if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
	  if (abs(rslo)>1.e-5) then
	    ff=tan(rphi)/tan(rslo)
	  else
!  set factor of safety to fixed value for flat slopes 	  
	    ff=finf
	  end if
	  zns=float(nzs)
	  zinc=(zmax(i)-zmin)/zns
	  z=zmin
	  fmn=1.e25
	  rf=0.0
	  z_loop: do j=1,nzs+1
	    znew=z
	    if(znew < 1.0e-30) znew =1.0e-30
	    if (abs(a1)>1.e-5) then
	      fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
	    else
	      fc(j)=0.d0
	    end if
	    pzero(j)=beta*(z-depth(i))
	    if (abs(z)>0.) then ! Formulas in next 2 lines apply only for z>0, test added 21 Feb 2013, RLB
	      zstar=z**2/dhat ! formula simplified 8 Jan 2013, RLB
	      tstar=t/zstar
	    end if
	    rf=0.0
	    temporal_loop: do n=1,nper
	      if(z==0.) then ! exact formula added for case of z=0, 8 Jan 2013, RLB
	        tdif1=t-capt(n)
	        if(tdif1>0.) then
	          rfa=sqrt(tdif1*dhat/pi)
	        else
	          rfa=0.0
	        end if
	        tdif2=t-capt(n+1)
	        if(tdif2>0.) then
	          rfb=sqrt(tdif2*dhat/pi)
	        else
	          rfb=0.0
	        end if
	      else ! z>0
	        captstar1=capt(n)/zstar
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
	        captstar2=capt(n+1)/zstar
	        tdif2=tstar-captstar2
	        if(tdif2 > 0.0) then	
	          x3=1./tdif2
	          x4=1./(sqrt(tdif2))
	          rf3=sqrt(1./(x3*pi))*exp(-x3)
	          rf4=derfc(x4)
	          rfb=rf3-rf4
	        else
	          rfb=0.0
 	        end if
 	      end if
	      rf=rf+rik(i+(n-1)*imax)*(rfa-rfb)
   	      if(rfa==0.0 .and. rfb==0.0) exit ! skip unnecessary cycles, RLB, 2/19/2015
	    end do temporal_loop
  	    ptran(j)=z*rf
            if(z==0) ptran(j)=rf ! formula for ptran at z=0 not normalized, 9 Jan 2013, RLB 
 	    p(j)=pzero(j)+ptran(j)
	    bline(j)=z*beta
	    ptest=p(j)-bline(j)
	    if(ptest > 0.0) then
	      p(j)=bline(j)
	    end if
	    if (abs(a1)>1.e-5) then
	      if(lpge0 .and. p(j)<0.) then !option added 4/15/2010
	        fw(j)=0.d0
	      else if (z>0.) then ! Added z>0 condition 2/12/2013, RLB
	        fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
	      end if
	    else
	      fw(j)=0.d0
	    end if
	    z=z+zinc
	  end do z_loop
! find new height of rising water table in zones of upward seepage   
	  if(rikzero(i)<0.0) then
 	    zinc=(zmax(i)-zmin)/zns
	    z=zmin
	    newdep=0.0
   	    z_loop_a: do j=1,nzs+1
	      if(p(j)<0.0) newdep=z
	      z=z+zinc
   	    end do z_loop_a
! adjust presures 
   	    z=zmin
	    z_loop_b: do j=1,nzs+1
	      if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
	      if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
	      z=z+zinc	    
	    end do z_loop_b
   	    end if
   	    z=zmin
	  fs_loop: do j=1,nzs+1
	    fs=ff+fw(j)+fc(j)
! frictional strength cannot be less than zero 	
	    if ((ff+fw(j))<0.) fs=fc(j)
 	    if (fs>finf) fs=finf
 	    if (z<=1.e-02) fs=finf 
!  	    if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
!  	    if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),pzero(j),ptran(j),&
!             & bline(j),fs
!  	    if (flag==-3) write(u1,'(6(g12.5,1x):)') z,p(j),fs ! added 4/14/2010 RLB, Revised 12/23/2010
	    if (fs<fmn) then
	      fmn=fs
	      zfmin(i)=z
	      pmn=p(j) ! revised 4/15/2010
	    end if
! Store pressure head and related output in 3-d arrays. Added 17Nov2014, RLB
 	        if(flag<0 .or. outp(1)) then
                  p3d(i,j)=p(j)
                  newdep3d(i)=newdep
                  dh3d(i)=0.d0
                end if 
 	    if(flag==-1) fs3d(i,j)=fs
 	    if(flag==-2) then
 	      fs3d(i,j)=fs
 	      ptran3d(i,j)=ptran(j)
 	      pzero3d(i,j)=pzero(j)
            end if
 	    if(flag==-3) then
 	      fs3d(i,j)=fs
! 	      th3d(i,j)=ths(zo(i))
            end if
            if(flag<=-4 .or. outp(1)) th3d(i,j)=ths(zo(i))
	    z=z+zinc
	  end do fs_loop
!  next statement assumes that computations begin at surface and work downward   
	  fsmin(i)=fmn
	    if(fmn==finf) then ! Added 30 Jan 2013, RLB 
	      pmn=p(nzs+1)
	      zfmin(i)=zmax(i)
            end if 
	  if(lpge0 .and. pmn<0.) then !option added 4/15/2010
	    pmin(i)=0.
	  else
	    pmin(i)=pmn
	  end if
!	  if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
!  	      dcf=0. ! this line used only for saturated infiltration model
!  	      chi=1.d0 !added 12/23/2010
!	      if (flag>=-6) call svijz(i,1,0.d0,newdep,ulog)
!	      if (flag<=-7 .and. flag>=-9) call svxmdv(i,1,0.d0,newdep,ulog)  ! Added 2/10/2012
!	  end if
	  if (mod(i,2000)==0) write (*,fmt='(2x,i10,2x,a1)',advance="no") i,char(13)
   	end do grid_loop
   	write(*,*)
   	write(*,*) imax, ' cells completed' 
   	write(ulog,*) imax, ' cells completed'
	return
	end subroutine iverson
