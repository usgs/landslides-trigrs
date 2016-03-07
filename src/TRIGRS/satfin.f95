	subroutine satfin(imx1,ulog,u1,nccs)
! 7/19/2006 Rex L. Baum, USGS 
!  computations for time series in the fully saturated case.	
	use grids; use input_vars 
	use model_vars
	use input_file_defs
	implicit none
	integer::i,j,jf,k,ulog,u1,imx1,nccs
	integer::nmn1,nmin1,nmax3,nmxs,nmns,nmax0
	logical:: lcvs 
	real:: qbij(nts+1)
	real (double)::rf(nzs+1),finf
	nmax3=0;nmax0=0;nmxs=0 ! modified 10 Nov 2014, RLB
	nmn1=nmax+1;nmin1=nmax+1
	write(ulog,*) 'Starting computations'
	write(ulog,*) 'for finite-depth saturated zone'
	write(*,*) 'Starting computations'
	write(*,*) 'for finite-depth saturated zone'
	write(*,*) 'Cells completed: '
! loop over all grid cells
	finf=10.
	  do i=1,imx1 
	    if (mod(i-1,2000)==0) write (*,fmt='(2x,i10,2x,a1)',advance="no") i-1,char(13) ! cells completed
	    if(slo(i)<slomin .or. slo(i)>slomax .or. zmax(i)<=0.0001) then ! default values for gently or steeply sloping cells 
	      do jf=1,nout
	        fsmin(i+(jf-1)*imax)=finf+1.
	        zfmin(i+(jf-1)*imax)=zmax(i)
	        pmin(i+(jf-1)*imax)=0.
	      end do
	      cycle
	    end if
	    lcvs=.true.
	    q=0.
	    do j=1,kper
	      if(j>nper) then
	        q(j)=0.
	      else
	        q(j)=ks(zo(i))*rik(i+(j-1)*imax)
	        if(q(j)>ks(zo(i))) write (ulog,*) '*q>Ks!', i,j,q(j),ks(zo(i))
	      end if
	    end do
!  use surface flux in infiltration computations
	    qb=0. ! initialize qb for case where ts>capt(nper+1)
	    ts=0.
	    do j=1,nts+1
    	      do k=1,kper
    	        if(ts>=capt(k) .and. ts<=capt(k+1)) qb(j)=q(k)
              end do
   	      if(outp(7)) rik1(i+(j-1)*imax)=qb(j)/ks(zo(i))
	      tcap(j)=ts ! pass to diffusion subroutine
   	      ts=ts+tinc_sat(j) ! Revised 25 June 2013 RLB
   	    end do
   	    do j=1,nts+1
	      qbij(j)=qb(j)/ks(zo(i))
  	    end do
 	        rf=0. 
	        call svgstp(u1,qbij,ulog,i,rf,nccs,lcvs,nmxs)
                nmns=nmn
	  end do
	  write(*,*)
	  write (*,*) imx1, ' cells completed' 
	  write (ulog,*) imx1, ' cells completed' 
!
	if(nmns>nmxs) nmns=nmxs
	write(ulog,*) 'Convergence data for saturated zone:'
	write(ulog,*) 'Max. terms used by sat-fin infinite series', nmxs
	write(ulog,*) 'Min. terms used by sat-fin infinite series', nmns
	write(ulog,*) 'Saturated-zone nonconvergent cells: '
	write(ulog,*) nccs
	return
	end subroutine satfin
