!  finite depth diffusion solution for rain infiltration.
!  by W.Z. Savage, Spring 2001, with modifications by R.L. Baum
!  (both) USGS, Latest revison, 4-6 Jan. 2016, RLB, corrected treatment of time-independent term in later-time formula
	subroutine savage(imx1,u1,profil,ulog,nccs)
	use grids; use input_vars
	use model_vars
	implicit none
	integer :: j,i,imx1,u1,nmx,ulog,n1,nccs,n,m 
	real (double):: finf,t1,t2,term1,term2,rn,znew
	real (double):: fierfc1,fierfc2,fierfc3,fierfc4
	real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
	real (double):: ar1,ar2,ar3,ar4,fs
	real (double):: rfa,rfb,rf,derfc,ff,rslo,rphi,fmn,ptest,pmn ! pmn added 4/15/2010
	real (double):: tol,delt1,delt2,t1old,t2old,tfac1,dif1 ! ,tfac2
	real (double):: ddg2rad,dlz,flt1,flt2,temp1,temp2,temp3,late_t,ksat
	real (double):: riksum ! added 5 Jan 2016, RLB
	logical :: lcv 
	character (len=255)::  profil
	write(ulog,*) 'Starting saturated-zone'
	write(ulog,*) 'computations for finite-depth'
	write(*,*) 'Starting saturated-zone'
	write(*,*) 'computations for finite-depth'
	pi=3.141592653589793
	ddg2rad=pi/180.D0
!  maximum value of Factor of Safety
	finf=10.
	late_t=5.0 ! optimized by trial and error 02/23/2015, RLB
!  loop steps through all grid cells
	write(*,*) 'Cells completed: '
	nmx=0
!	nmn=1+mmax ! initalization of nmn moved to main program 07 Jan 2013, RLB 
	grid_loop: do i=1,imx1
	  rslo=slo(i)
	  if(rslo<slomin .or. rslo>slomax .or. zmax(i)<=0.0001) then
	    fsmin(i)=finf+1.
	    zfmin(i)=zmax(i)
	    pmin(i)=0.
	    if (mod(i,2000)==0) write (*,fmt='(2x,i10,2x,a1)',advance="no") i,char(13)
	    cycle grid_loop
	  end if
	  rphi=phi(zo(i))
	  a1=sin(rslo)
	  b1=cos(rslo)
	  dif1=dif(zo(i))/(b1*b1)
	  ksat=ks(zo(i))
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
	  else
!  set factor of safety to fixed value for flat slopes 	  
	    ff=finf
	  end if
	  zns=float(nzs)
	  zinc=(zmax(i)-zmin)/zns
	  z=zmin
	  lcv=.true.
! compute depth of saturated zone, dlz
	  dlz=zmax(i) 
	  Z_loop: do j=1,nzs+1
	    znew=z
	    if(znew < 1.0e-30) znew =1.0e-30
	    if (abs(a1)>1.e-5) then
	      fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
	    else
	      fc(j)=0.d0
	    end if
	    pzero(j)=beta*(z-depth(i))
	    rf=0.0
            riksum=rik(i);  ! added 5 Jan 2016, RLB
            temp2=dlz*(3.*(dlz-z)*(dlz-z)-dlz*dlz)/(6.*dlz*dlz) !time-independent factor, Moved 5 Jan 2016, RLB  
	    temporal_loop: do m=1,nper
	      tdif1=t-capt(m)
	      tfac1=tdif1*dif1/(dlz*dlz) ! Moved to top level of temporal_loop from line 89, 4 Jan 2016, RLB
! the next if block assumes that no change in transient pore pressure
! occurs if no flux at the top of the saturated zone.
!	      if(rik(i+(m-1)*imax)==0.) then
!	        cycle temporal_loop	      
!	      end if
	      if(tdif1 > 0.0) then
! corrected diffusivity term in next line (divide by b1*b1) 	      
                t1=sqrt(tdif1*dif(zo(i))/(b1*b1))
	        if (t1<1.0e-29) t1=1.0e-29
	        term1=0.0
	          if(tfac1>late_t) then !Added later time solution to improve performance and eliminate nonconvergent cells, 2/23-24/2015, RLB 
	            temp1=tdif1*dif1/dlz ! Corrected 7 Jan 2016
                            if(m>1) riksum=riksum+rik(i+(m-1)*imax)-rik(i+(m-2)*imax) ! Added 8 Jan 2016, RLB
	            rf=rf+rik(i+(m-1)*imax)*temp1+riksum*temp2 ! revised 8 Jan 2016
	            series_a_lt: do n=1,mmax ! later time solution
	              rn=float(n)
	              ar1=rn*pi*(dlz-z)/(dlz)
	              ar2=(rn*pi)/dlz
	              flt1=exp(-dif1*tdif1*ar2*ar2)*cos(ar1)
	              flt1=float(-1**n)/(rn*rn)*flt1
! test for convergence of series to within 1/1e+06 of previous value	      
	              t1old=term1
	              tol=abs(term1/1e+06) ! tol=term1/10000. revised 2/17/2015, RLB
   	              term1=term1+flt1
 	              delt1=abs(term1-t1old)
 	              n1=n
 	              if(delt1<=tol) exit
   	            end do series_a_lt
	            rfa=2.*dlz*term1/(pi*pi) 
                  else
	            series_a_et: do n=1,mmax ! early time solution
	              rn=float(n)
	              ar1=((2.*rn-1.)*dlz-(dlz-z))/(2.*t1)
	              ar2=((2.*rn-1.)*dlz+(dlz-z))/(2.*t1)
	              fierfc1=exp(-ar1**2)/sqrt(pi)-ar1*derfc(ar1)
	              fierfc2=exp(-ar2**2)/sqrt(pi)-ar2*derfc(ar2)
! test for convergence of series to within 1/1e+06 of previous value	      
	              t1old=term1
	              tol=abs(term1/1e+06) ! tol=term1/10000. revised 2/17/2015, RLB
   	              term1=term1+fierfc1+fierfc2
 	              delt1=abs(term1-t1old)
 	              n1=n
 	              if(delt1<=tol) exit
   	            end do series_a_et
	            rfa=2.*t1*term1 
                  end if
                  if(lcv .and. delt1>tol) then
                    nccs=nccs+1
                    nv(i)=1
                    lcv=.false.
                  end if
  	          if(n1>nmx) nmx=n1
   	          if(n1<nmn) nmn=n1
	        else
	          rfa=0.0
	        end if
	        tdif2=t-capt(m+1)
!                        tfac2=tdif2*dif1/(dlz*dlz) ! Moved up from line 142, 4 Jan 2016, RLB
	        if(tdif2 > 0.0) then
!  corrected diffusivity term in next line (divide by b1*b1) 	      
                  t2=sqrt(tdif2*dif(zo(i))/(b1*b1))
	          if (t2<1.0e-29) t2=1.0e-29
	          term2=0.0 
	          if(tfac1>late_t) then !
	            temp3=tdif2*dif1/dlz	          
	            rf=rf-rik(i+(m-1)*imax)*(temp3)
	            series_b_lt: do n=1,mmax ! later time solution
	              rn=float(n)
	              ar3=rn*pi*(dlz-z)/(dlz)
	              ar4=(rn*pi)/dlz
	              flt2=exp(-dif1*tdif2*ar4*ar4)*cos(ar3)
	              flt2=float(-1**n)/(rn*rn)*flt2
! test for convergence of series to within 1/1e+06 of previous value	      
	              t2old=term2
	              tol=abs(term2/1e+06) ! tol=term2/10000. revised 2/17/2015, RLB
   	              term2=term2+flt2
 	              delt2=abs(term2-t2old)
 	              n1=n
 	              if(delt2<=tol) exit
   	            end do series_b_lt
	            rfb=2.*dlz*term2/(pi*pi) 
                  else
	            series_b_et: do n=1,mmax ! early time solution
	              rn=float(n)
	              ar3=((2.*rn-1.)*dlz-(dlz-z))/(2.*t2)
	              ar4=((2.*rn-1.)*dlz+(dlz-z))/(2.*t2)
	              fierfc3=exp(-ar3**2)/sqrt(pi)-ar3*derfc(ar3)
	              fierfc4=exp(-ar4**2)/sqrt(pi)-ar4*derfc(ar4)
! test for convergence of series to within 1/1e+06 of previous value	      
	              t2old=term2
	              tol=abs(term2/1e+06) ! tol=term2/10000. revised 2/17/2015, RLB
   	              term2=term2+fierfc3+fierfc4
 	              delt2=abs(term2-t2old)
 	              n1=n
 	              if(delt2<=tol) exit
   	            end do series_b_et
	            rfb=2.*t2*term2 
                  end if
   	          if(lcv .and. delt2>tol) then
   	            nccs=nccs+1
                    nv(i)=1
   	            lcv=.false.
   	          end if
   	          if(n1>nmx) nmx=n1
   	          if(n1<nmn) nmn=n1
	        else
	          rfb=0.0
	        end if
	        if(tfac1>late_t) then
   	          rf=rf-rik(i+(m-1)*imax)*(rfa-rfb)
                else
   	          rf=rf+rik(i+(m-1)*imax)*(rfa-rfb)
                end if
   	      if(rfa==0.0 .and. rfb==0.0) exit ! skip unnecessary cycles, RLB, 2/19/2015
   	    end do temporal_loop 
            ptran(j)=rf
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
   	  z=zmin
	  fmn=1.e25
	  Z_FS_loop: do j=1,nzs+1
	    fs=ff+fw(j)+fc(j)
! frictional strength cannot be less than zero	
	    if ((ff+fw(j))<0.) fs=fc(j)
 	    if (fs>finf) fs=finf
 	    if (z<=1.e-02) fs=finf 
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
   	  end do Z_FS_loop
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
	  if (mod(i,2000)==0) write (*,fmt='(2x,i10,2x,a1)',advance="no") i,char(13)
   	end do grid_loop
   	write(*,*)
   	write(*,*) imx1, ' cells completed'
   	write(ulog,*) imx1, ' cells completed'
	if(t==0 .and. nper==1) then ! added 4/14/2010
	  nmx=0; nmn=0
	end if
	write(ulog,*) 'Max. terms used by sat-fin infinite series', nmx
	write(ulog,*) 'Min. terms used by sat-fin infinite series', nmn
	write(ulog,*) 'Saturated-zone nonconvergent cells: '
	write(ulog,*) nccs
	return
	end subroutine savage
