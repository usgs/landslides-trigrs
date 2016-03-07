	subroutine unsinf(imx1,ulog,u1,ncc,nccs)
! By Rex L. Baum, USGS, Latest revision, 27 Jun 2013 
	use grids; use input_vars 
	use model_vars
	use input_file_defs
	implicit none
	integer::i,j,jf,k,ulog,u1,imx1,ncc,nccs,nmax0 ! ,nmx
	integer::nmn1,nmin1,nmax3 !,nmns,nmxp,nmnp,nmxs
	integer:: ivctr
	logical:: lcv,lwt ! ,lcvs
	real:: delwt,dwt,zwt,qbij(nts+1) 
        real::testqk,tolqk ! Added 2/2/2011 RLB
	real (double)::rf(nzs+1),finf,vqt,qta,al,qzmax 
	real (double)::ddwt,sqin,intq(nts+1),b,dhwt(nts+1),delh 
	real (double)::qtn(2*nts+1),intq1(nts+1),vqtn,cd
	nmax3=0;nmax0=0
	nmn1=nmax+1;nmin1=nmax+1
        ivctr=0 !nmxp=0;nmnp=0;   ! Added 29 Jan 2013, RLB 
	intq=0.d0; intq1=0.d0 ! Added 29 Jan 2013, RLB 
	write(ulog,*) 'Starting coupled saturated & unsaturated zone'
	write(ulog,*) 'computations for infinite-depth saturated zone'
	write(*,*) 'Starting coupled saturated & unsaturated zone'
	write(*,*) 'computations for infinite-depth saturated zone'
	write(*,*) 'Cells completed: '
! loop over all grid cells
	finf=10.
	  grid_loop: do i=1,imx1 
	    if (mod(i-1,2000)==0) write (*,fmt='(2x,i10,2x,a1)',advance="no") i-1,char(13) ! cells completed
	    if(slo(i)<slomin .or. slo(i)>slomax .or. zmax(i)<=0.0001) then ! default values for gently or steeply sloping cells 
	      do jf=1,nout
	        fsmin(i+(jf-1)*imax)=finf+1.
	        zfmin(i+(jf-1)*imax)=zmax(i)
	        pmin(i+(jf-1)*imax)=0.
	      end do
	      cycle
	    end if
	    lcv=.true. !;lcvs=.true.
	    q=0.;qb=0 ! qb initialization added 29 Jan 2013, RLB 
	    tolqk=ks(zo(i))*5.e-07 ! Moved 29 Jan 2013, RLB 
	    do j=1,kper
	      if(j>nper) then
	        q(j)=0.
	      else
! 8/18/2009 RLB added optional offset of background flux to prevent excessive drying during periods of zero infiltration.	      
	        if(bkgrof) then
	          q(j)=ks(zo(i))*(rik(i+(j-1)*imax)+rikzero(i))
! RLB 2/2/2001 revised test
	          testqk=q(j)-(ks(zo(i))+rizero(i))
!	          tolqk=ks(zo(i))*5.e-07
	          if(testqk>tolqk) write (ulog,*) '*q>Ks+ri!', i,j,q(j),ks(zo(i))+rizero(i) ! Added 2/2/2011 RLB
	        else
	          q(j)=ks(zo(i))*rik(i+(j-1)*imax)
	          testqk=q(j)-ks(zo(i))
	          tolqk=ks(zo(i))*5.e-07          
	          if(testqk>tolqk) write (ulog,*) '*q>Ks!', i,j,q(j),ks(zo(i)) ! Moved 2/2/2011 RLB
	        end if
	      end if
	    end do
	    qmax=maxval(q)
	    b=cos(slo(i))
! next lines compute depth to top of capilarly fringe.
	    if(unsat(zo(i))) then
	      dcf=depth(i)-1.d0/(alp(zo(i)))
	    else
	      dcf=0.
	    end if 
	    if(lps0 .and. unsat(zo(i))) then
	      dusz=depth(i)-1.d0/(alp(zo(i)))
	    else 
	      dusz=depth(i)
	    end if
! set value of beta (Iverson's beta line)
! and maximum drainage rate at water table, qzmax	    
	    cd=1.0d0 ! complete drainage for infinite depth basal flow boundary	         
	    select case (flowdir)
	      case ('slope')
	      beta=b*b
	      qzmax=(1.d0-beta)*cd*ks(zo(i)) 
	      case ('hydro')
	      beta=1.d0
	      qzmax=0.d0 
	      case default
	      beta=b*b-rikzero(i) ! 2/12/09 corrected formula for beta & qzmax
	      qzmax=(1.d0-beta)*cd*ks(zo(i))-cd*rizero(i) 
	    end select
	    delwt=0.;dwt=depth(i);zwt=depth(i);ddwt=depth(i)
	    lwt=.false.
	    ts=tmin
	    if(dcf>0. .and. (unsat(zo(i)) .or. igcap(zo(i)))) then ! updated to enforce non-zero depth, 4/18/2013, RLB
! compute flux and pore-pressure rise at each time step	   
	      vqt=0.;vqtn=0.;qta=rizero(i);sqin=0.
 	      al=alp(zo(i))*b*b
	      if(outp(8)) write (ulog,*) 'ts,    qt '! times and basal flux to log file
	      call roots(nmax,r,dusz,al,eps,pi)
	      flux_loop: do j=1,2*nts+1
	        call flux(i,kper,ts,j,lcv,ncc,nvu(i),lwt) 
                 if(nmax1>nmax2) nmax2=nmax1
                 if(nmn<nmin) nmin=nmn
                 if(qt<rizero(i)) qt=rizero(i)
! RLB 2/3/2011, Added case for bkgrof=.true.               
	         if(bkgrof) then
	           if(qt>ks(zo(i))+rizero(i)+tolqk) then
	             write(ulog,*) 'Error! Basal flux exceeds Ks at'
	             write(ulog,*) 'cell ',i, ', timestep ',j 
	             write(ulog,*) 'flux ',qt, ', Ks ',ks(zo(i))+rizero(i) 
	             qt=ks(zo(i))
	           end if
	         else
	           if(qt>ks(zo(i))+tolqk) then
	             write(ulog,*) 'Error! Basal flux exceeds Ks at'
	             write(ulog,*) 'cell ',i, ', timestep ',j 
	             write(ulog,*) 'flux ',qt, ', Ks ',ks(zo(i))
	             qt=ks(zo(i))
	           end if
	         end if
	         if(outp(8)) write (ulog,*) ts,qt ! times and basal flux to log file
! drain off excess basal flux
	         if(qt>qzmax) then
	           qtime(j)=qt-qzmax
	           else
	           qtime(j)=0.d0
	         end if
	        qtn(j)=qt
	        qta=qt
   	        ts=ts+tinc/2.d0
	      end do flux_loop
	      call dsimps(nts,tinc/2.d0,qtime,intq)
	      call dsimps(nts,tinc/2.d0,qtn,intq1)
	      if(outp(8)) write (ulog,*) 'Time, Cumulative volume in, Cumulative background flux,&
              & Cumul. volume out,Cuml. absorbed, Cuml. qin-qout, qout not drained, Water table rise'
	      ts=tmin
	      wt_rise_loop: do j=1,nts+1
	        jf=jsav(j)
	        rf=0.0 
	        vqt=intq(j)-ts*rizero(i)
	        if(vqt<0.) vqt=0.d0
	        vqtn=intq1(j)-ts*rizero(i)
	        sqin=0.
	        sum_q_in_loop: do k=1,nper
	          if(ts>capt(k) .and. ts<=capt(k+1)) then
    	            qts(j)=q(k) !; write(*,*) 'j,k,qts(j),q(k),ts,capt(k) ', j,k,qts(j),q(k),ts,capt(k) ! Added 21 Feb 2013, RLB 
	          endif
	          if(ts>=capt(k+1)) then
	            sqin=sqin+(capt(k+1)-capt(k))*q(k)
	          end if
	          if(ts>capt(k) .and. ts<capt(k+1)) then
	            sqin=sqin+(ts-capt(k))*q(k)
	          end if
	        end do sum_q_in_loop
	        if(jf>0 .and. lskip) then
! compute unsaturated zone pressure & water content
	          call unsth(i,j,ncc,kper,ts,nmax0,& 
	          &lcv,ulog,vqt,delh,nmn1,sqin,vqtn)
                  dhwt(j)=delh
                  if(nmax0>nmax3) nmax3=nmax0
                  if(nmn1<nmin1) nmin1=nmn1
                else if(lskip) then
                  continue
                else
	          call unsth(i,j,ncc,kper,ts,nmax0,& 
	          &lcv,ulog,vqt,delh,nmn1,sqin,vqtn)
                  dhwt(j)=delh
                  if(nmax0>nmax3) nmax3=nmax0
                  if(nmn1<nmin1) nmin1=nmn1
	        end if	        
	          tcap(j)=ts ! pass to diffusion subroutine
	          tcap(j+1)=ts+tinc
   	          ts=ts+tinc
	        if(jf>0 .and. lskip) then
! compute pressure diffusion in saturated zone
	          call pstpi(u1,dhwt,dwt,ulog,i,j-1,rf,tcap(j),jf) ! Revised 29 Jan 2013, RLB 
!                  nmxp=nmx;nmnp=nmn
	          if(j>1) then
! Check change in water table depth and adjust dusz if needed
	            delwt=abs(dwt-zwt)*1000.
	            if(delwt>dwt) then
	              lwt=.true.
	              dwt=zwt
	            end if
	          end if   
	        else if(lskip) then
	          continue
	        else
! compute pressure diffusion in saturated zone
	          call pstpi(u1,dhwt,dwt,ulog,i,j-1,rf,tcap(j),jf) ! Revised 29 Jan 2013, RLB 
!                  nmxp=nmx;nmnp=nmn
	          if(j>1) then
! Check change in water table depth and adjust if needed
	            delwt=abs(dwt-zwt)*1000.
	            if(delwt>dwt) then
	              lwt=.true.
	              dwt=zwt
	            end if
	          end if   
	        end if	        
	      end do wt_rise_loop
! map unsaturated zone outflux to grid
! there are 2*nts+1 increments in qtime()
	      do k=1,nts  
   	        qb(k)=qtime(2*k+1)
   	        if(outp(7)) rik1(i+(k-1)*imax)=qb(k)/ks(zo(i))
   	      end do
	    else ! top of capillary fringe at ground surface, so use surface flux
	      qb=0. ! initialize qb for case where ts>capt(nper+1)
	      dwt=depth(i)
	      delh=0.;rf=0.;zwt=depth(i)
	      do j=1,nts+1
    	        do k=1,kper
    	          if(ts>=capt(k) .and. ts<=capt(k+1)) qb(j)=q(k)
                end do
   	        if(outp(7)) rik1(i+(j-1)*imax)=qb(j)/ks(zo(i))
	        tcap(j)=ts ! pass to diffusion subroutine
   	        ts=ts+tinc
   	      end do
   	      do j=1,nts+1
	        qbij(j)=qb(j)/ks(zo(i))
  	      end do
	          rf=0.
	          ivctr=ivctr+1
 	          call ivestp(u1,qbij,ulog,i,rf) ! Revised 29 Jan 2013, RLB 
!                  nmns=nmn
	    end if
	  end do grid_loop
	  write(*,*)
	  write (*,*) imx1, ' cells completed' 
	  write (ulog,*) imx1, ' cells completed' 
!
	if(nmin>nmax2) nmin=nmax2; if(nmin1>nmax3) nmin1=nmax3
	write(ulog,*) 'Convergence data for unsaturated zone:'
	write(ulog,*) 'Maximum terms used by Fourier series', nmax2,nmax3
	write(ulog,*) 'Minimum terms used by Fourier series', nmin,nmin1
	write(ulog,*) 'Unsaturated zone nonconvergent cells: '
	write(ulog,*) ncc
!	if(nmnp>nmxp) nmnp=nmxp; if(nmns>nmxs) nmns=nmxs
	write(ulog,*) 'Convergence data for saturated zone:'
! ivestp() and pstpi() use single-term solutions to compute pressure head--Convergence is not an issue 	
	write(ulog,*) 'Terms used by infinite-depth '
	write(ulog,*) 'solutions ivestp() and pstpi():', 1, 1
	write(ulog,*) 'Saturated-zone nonconvergent cells: '
	write(ulog,*) nccs !initialized to zero in trigrs main
	write(ulog,*) 'Cells using ivestp() = ',ivctr
	return
	end subroutine unsinf
