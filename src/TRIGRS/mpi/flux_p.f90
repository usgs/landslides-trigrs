! by W.Z. Savage with modifications by R.L. Baum, Latest revision 29 Jan 2013, RLB
!
! MPI calls & MPI-related changes by M. Alvioli, November 2014
!
subroutine flux_p(ic,iper,t0,j1,lcv,ncc,nv0,lwt)
  use grids; use input_vars
  use model_vars
  use mpi
  implicit none
  integer:: ic
  integer:: j1,iper,i,k,k1,ncc,nv0 !,i1 removed 12/28/2010
  real (double):: al,qa,tol
  real (double):: t0,tf,q1a,q1b,q2a,q2b,qta,qtb,qtop !,bot
  real (double):: q2old,delq2,tdif1,tdif2,tstar
  real (double):: q0a,q0b 
  real (double):: z,psih,qlb,ck,th,b1
  real (double):: bot(nmax),qtops(nmax),dseep ! added qtops & made arrays 12/28/2010
  !
  logical:: lcv,lwt
  integer myrank,isize,ierr
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  tol=1.0e-06 ! Revised 3/4/2015, RLB
  b1=cos(slo(ic))
  ! coordinate transformation, corrected 9/5/06	
  al=alp(zo(ic))*b1*b1
  qa=rizero(ic)
  tf=al*ks(zo(ic))/(ths(zo(ic))-thr(zo(ic)))
  qt=0.0; qta=0.0; qtb=0.0 
  nmax1=0
  nmn=nmax+1
  do k=1,nmax
     bot(k)=1.+al*dusz/2.+2.*al*dusz*r(k)**2
     qtops(k)=r(k)*sin(r(k)*al*dusz)
  end do
  flux_time_loop: do  i=1,iper
     !	  i1=1+(i-1)*tx
     tdif1=t0-capt(i)
     if (tdif1 < 0.) exit ! jump out of loop rather than compute extra zeros.
     if(tdif1 > 0.0) then
        tstar=tf*tdif1
        if(tstar < smt) then
           ! Early-time solution (ETS) ...............
           z=dusz
           qlb=ks(zo(ic))
           ! Corrected ths, thr to include index 12/20/2010, RLB		
           call smallt(tstar,al,qa,q(i),dusz,qlb,ths(zo(ic)),thr(zo(ic)),qta,psih,ck,th,z)
           nmn=1
           if(outp(8).and.myrank.eq.0) write(*,*) 'cell, time, t* ', ic,t0,tstar, &
                &' Using ETS for basal flux' !Revised 2/2/2011 RLB
        else    
           ! later-time solution ...............	    
           q2a=0.0
           do k=1,nmax
              qtop=qtops(k)*exp(-(r(k)**2)*tstar) !12/28/2010
              q2old=q2a
              q2a=q2a+qtop/bot(k)  !12/28/2010
              delq2=abs(q2a-q2old)/ks(zo(ic)) ! Revised 23 Oct 2013, RLB
              k1=k
              if(abs(q2a)<=tol .and. k>3) exit 
              if(delq2<=tol) exit  ! Revised 3 Mar 2015, RLB
           end do
           if(lcv) then
              if((delq2>tol) .and. (abs(q2a)>tol .and. k>3)) then  ! Revised 3 Mar 2015, RLB
                 ncc=ncc+1
                 nv0=1
                 lcv=.false.
              end if
           end if
           if(k1>nmax1) nmax1=k1
           if(k1<nmn) nmn=k1
           q0a=q(i)
           q1a=4.d0*(q(i)-qa)*exp(al*dusz/2.)
           qta=q0a-q1a*q2a*exp(-tstar/4.)
        end if
     else 
        qta=0.0
     end if
     tdif2=t0-capt(i+1)
     if(tdif2 > 0.0) then
        tstar=tf*tdif2
        if(tstar < smt) then
           ! Early-time solution ...............	    
           z=dusz  ! inititalizations added 29 Jan 2013, RLB 
           qlb=ks(zo(ic))
           ! Corrected ths, thr to include index 12/20/2010, RLB		
           call smallt(tstar,al,qa,q(i),dusz,qlb,ths(zo(ic)),thr(zo(ic)),qtb,psih,ck,th,z)
           nmn=1
           if(outp(8).and.myrank.eq.0) write(*,*) 'cell, time, t* ', ic,t0,tstar, &
                &' Using ETS for basal flux' !Revised 2/2/2011 RLB
        else
           ! later-time solution ...............	    
           q2b=0.0
           do  k=1,nmax
              qtop=qtops(k)*exp(-(r(k)**2)*tstar) !12/28/2010
              q2old=q2b
              q2b=q2b+qtop/bot(k) !12/28/2010
              !	        delq2=abs((q2b-q2old)/q2b)
              delq2=abs(q2b-q2old)/ks(zo(ic))  ! Revised 23 Oct 2013, RLB
              k1=k
              if(abs(q2b)<=tol .and. k>3) exit 
              if(delq2<=tol) exit  ! Revised 3 Mar 2015, RLB
           end do
           if(lcv) then
              if((delq2>tol) .and. (abs(q2b)>tol .and. k>3)) then  ! Revised 3 Mar 2015, RLB
                 ncc=ncc+1
                 nv0=1
                 lcv=.false.
              end if
           end if
           if(k1>nmax1) nmax1=k1
           if(k1<nmn) nmn=k1
           q0b=q(i)
           q1b=4.d0*(q(i)-qa)*exp(al*dusz/2.)
           qtb=q0b-q1b*q2b*exp(-tstar/4.)
        end if
     else
        qtb=0.0
     end if
     qt=qt+qta-qtb
  end do flux_time_loop
  if(qmax<qt .or. qt<0) then ! moved after end of loop, 7 Jan 2013
     dseep=ks(zo(ic))/(ths(zo(ic))*t0) !Added check of distance traversed at saturated seepage velocity 12/28/2010 RLB
     if(dseep < dusz) then
        qt=0.d0
     else
        if(myrank.eq.0) then
           write(*,*) 'Error computing basal flux!' !Revised 12/28/2010, RLB
           write(*,*) 'Cell, Depth, Time, t*, Max. Input flux, Basal flux: ',ic,dusz,t0,tstar,qmax,qt
           write(*,*) ''
        endif
     end if
  end if
  return
end subroutine flux_p
