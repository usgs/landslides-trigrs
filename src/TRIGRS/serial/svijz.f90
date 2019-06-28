!  Routine to prepare pressure head & water content (thz) data for export to ijz format
!  Calling routine should set value of dcf if it is not initialized in main.
!  Rex L. Baum, 25 Feb 2010, latest revision 22 Aug 2013
! Volumetric water content substituted for relative saturation throughout 8/22/13, RLB
subroutine svijz(i,jf,delh,newdep,ulog)
use grids; use input_vars 
use model_vars
use input_file_defs
implicit none
integer::i,jf,ulog,m,m1,nw,mfrst,mlast!,mnp(1),mctr,mmid
integer::mrpx(1),mrpn(1),dinc,wtctr,wctr,wptr(nzs+1),wtptr(nzs+1)
real (double)::delh,newdep,zns,zinc,fdinc !,pab(nzs+1),pmn
real (double):: dwat,z,ptop,pbot,zwat(nzs+1),zbot,ztop,th_top,th_bot,zwat0
real (double):: zdeep,pdeep,th_deep,zrpn,zrpx,relhgt!,x,y,zmid,zmn(1),zmx(1)
real (double):: tol1,tol2,mp,plin(nzs+1),rp(nzs+1),rpmin,rpmax
logical ::ldeep, lopmx,lopmn,loxlo
!  compute water-table elevations from surface elevation and depths
zbot=elev(i)-zmax(i)
ztop=elev(i)-zmin
zwat0=elev(i)-depth(i) ! Initial water table depth
zwat=zwat0 ! initialize array values to initial water table depth
wtctr=0; wctr=0; wptr=0; wtptr=0
! Check for deepz <= 0 to skip deep node
ldeep=.true.
zdeep=elev(i)-deepz 
if(zdeep>ztop) ldeep=.false.
if(zdeep>zbot .and. ldeep) then
  zdeep=zbot-zmax(i) 
end if
zns=float(nzs) 
zinc=(zmax(i)-zmin)/zns
ptop=p(1)
pbot=p(nzs+1)
do nw=1, nzs+1 ! Adjust water content for grid cells where water table rise calculated with unsaturated model
  if(p(nw)>=0.d0) thz(nw)=ths(zo(i))
end do
th_top=thz(1) 
th_bot=thz(nzs+1)
th_deep=ths(zo(i))
! estimate height of water table from results of unsaturated infiltration model
if(dcf>0. .and. (unsat(zo(i)) .or. igcap(zo(i)))) then
  if(rikzero(i)>=0.0) then
    zwat(1)=elev(i)-(dusz-delh) ! statement in unsth() accounts for initial wt below zmax(i)
  else
    zwat(1)=elev(i)-newdep
  end if
wtctr=1;wtptr(wtctr)=1;wptr(wtctr)=1  ! changed wtptr(wctr) to wtptr(wtctr) 18 April 2013, RLB
else if (dcf<=0.) then
! find water table(s) from results of saturated infiltration model
  call dzero_brac(nzs,p,wtctr,wtptr) ! wtptr(j)=1 for intervals containing a water table
  thz=ths(zo(i))
  wctr=0;dwat=0.d0;wptr=0
  do nw=1,nzs+1
    if(wtptr(nw)>0) then
      if(wtptr(nw)==2) then ! water table at node
        wctr=wctr+1
        z=zmin+zinc*float(nw-1)
        dwat=z 
      end if
      if(wtptr(nw)==1) then ! water table between nodes
        wctr=wctr+1
        z=zmin+zinc*float(nw-1)
        dwat=z+zinc*(0.d0-p(nw))/(p(nw+1)-p(nw))
      end if
      if(rikzero(i)>=0.0) then
        zwat(nw)=elev(i)-dwat 
        wptr(wctr)=nw
      else
        zwat(nw)=elev(i)-newdep
        wptr(wctr)=nw
      end if
    end if
    if(wctr==wtctr) exit  ! all zero points accounted for
  end do
  if(p(nzs+1)<0.d0) then ! water table below bottom node (below zmax(i))
    wtptr(nzs+1)=1
    wtctr=wtctr+1;wctr=wctr+1
    zwat(nzs+1)=elev(i)-depth(i)+beta*(p(nzs+1)-p0zmx)  ! Approximate formula RLB 5/3/2013.
    if(zwat(nzs+1)<(elev(i)-depth(i))) zwat(nzs+1)=elev(i)-depth(i) ! prevents water table from falling below initial value
    wptr(wtctr)=nzs+1
  end if
end if
if(outp(1)) then
  if(wtctr==0)then
    if(p(1)>0. .and. zmin>=0.)then ! water table at ground surface and zmin is at or below surface.  Added 5/3/2013, RLB 
      wtctr=wtctr+1;wctr=wctr+1
      wptr(wtctr)=1
      zwat(1)=elev(i);p(1)=0.
    else 
      write(*,*) 'Error in svxmdv()!: wtctr not initialized at cell', i
      write(*,*) 'wptr(1), wtptr(1) ', wptr(1), wtptr(1)
    endif
  endif
  select case (el_or_dep) ! water table elevation or depth
    case('eleva')
      wtab(i+(jf-1)*imax)=zwat(wptr(wtctr)) ! use lowest computed water table
    case('depth')
      wtab(i+(jf-1)*imax)=elev(i)-zwat(wptr(wtctr))
    case default
      wtab(i+(jf-1)*imax)=zwat(wptr(wtctr)) ! Use lowest computed water table
  end select
end if
if(flag>-4 .or. flag<-7) return 
select case (deepwat) ! pore pressure at deep node for SCOOPS 
  case('zero') ! zero pressure
    pdeep=0.
  case('flow') ! pressure consistent with user specified flow direction, use lowest computed water table
    pdeep=beta*(elev(i)-depth(i)-zdeep) ! use initial depth to water table
  case('hydr') ! hydrostatic pressure option
    pdeep=(elev(i)-depth(i)-zdeep) ! use initial depth to water table
  case('relh') ! hydrostatic pressure reduced by relative height and constant factor
    relhgt=(elev(i)-zmn(1))/(zmx(1)-zmn(1))
    pdeep=(elev(i)-depth(i)-zdeep)*(1.-relhgt/3.) ! use initial depth to water table
  case default
    pdeep=0.
end select
! Save data in ijz format, full (flag==-4) or downsampled (flag==-5) depth profile
if(flag==-4 .or. flag==-5) then
  if(spcg<1) spcg=1
  dinc=1; if(flag==-5) dinc=spcg !dinc used for downsampled ijz file.  !Added spcg 2/14/2012 RLB
  fdinc=float(dinc)
  z=ztop
  if(dcf>0. .and. unsat(zo(i))) then !Unsaturated infiltration model
    do m=1,nzs+1,dinc
      write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),z,p(m),thz(m)
      if(zwat(1) < z .and. zwat(1) > (z-fdinc*zinc) .and. p(m)<0.) then ! eliminate repeated entries for water table
        if(zwat(1) > (zbot*1.00001)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & ix(i),jy(i),zwat(1),0.d0,ths(zo(i))
      end if
      z=z-(zinc*fdinc)
    end do
    if(zwat(1) < (zbot/1.00001)) then ! water table below zmax(i)
      write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zwat(1),0.d0,ths(zo(i)) 
    end if
    if(ldeep) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),zdeep,pdeep,th_deep
  else if (dcf<=0.) then ! Saturated infiltration model
    do m=1,nzs+1,dinc
      write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),z,p(m),thz(m)
      if(dinc==1) then
        if(wtptr(m)==1) then
          write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(m),0.d0,ths(zo(i))
        end if
      else if(dinc>1) then
        do m1=m, m+dinc-1
          if(wtptr(m1)==1) then
            write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & ix(i),jy(i),zwat(m1),0.d0,ths(zo(i))
          end if
        end do
      end if
      z=z-(zinc*fdinc)
    end do
    if(ldeep) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),zdeep,pdeep,th_deep
  end if
end if  
if(flag==-6) then ! Sparce output option
! find points at maximum (+) and minimum (-) excursions from linear trend between ground surface and basal water table
  lopmn=.false.; lopmx=.false.; loxlo=.false.
  plin=0.d0; rp=0.d0
  mp=ptop/(ztop-zwat(nzs+1))
  z=ztop
  plin(1)=ptop
  do m=2,nzs+1
    z=z-zinc
    if(z<zwat(1) .and. dcf>0. .and. unsat(zo(i))) exit
    plin(m)=ptop-mp*(ztop-z)
    rp(m)=p(m)-plin(m)
  end do
  rpmin=minval(rp); mrpn=minloc(rp)
  rpmax=maxval(rp); mrpx=maxloc(rp)
  tol1=abs(p(mrpn(1))/100.)
  tol2=abs(p(mrpx(1))/100.)
  zrpn=ztop-zinc*float(mrpn(1)-1)
  zrpx=ztop-zinc*float(mrpx(1)-1)
  if(abs(rpmin) >= tol1 .and. rpmin<0.d0 .and. zrpn>zwat(1)) lopmn=.true.
  if(abs(rpmax) >= tol2 .and. rpmax>0.d0 .and. zrpx>zwat(1)) lopmx=.true.
  mfrst=mrpn(1);mlast=mrpx(1)
  if(mrpn(1)>mrpx(1)) then
    loxlo=.true.
    mfrst=mrpx(1);mlast=mrpn(1)
  end if
! Save data in ijz format
  if(dcf>0. .and. unsat(zo(i))) then 
    write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),ztop,ptop,th_top !,'top1' 
    if(loxlo) then
      if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zrpx,p(mrpx(1)),thz(mrpx(1)) !,'max1'
      if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zrpn,p(mrpn(1)),thz(mrpn(1)) !,'min1'
    else
      if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zrpn,p(mrpn(1)),thz(mrpn(1)) !,'min2'
      if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zrpx,p(mrpx(1)),thz(mrpx(1)) !,'max2'
    end if
    if(zwat(1) < ztop .and. ptop < 0.) then ! eliminate repeated entries for water table Unsaturated model has only one water table
      if(zwat(1) > (zbot*1.00001)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zwat(1),0.d0,ths(zo(i)) !,'watu1'
    end if
    write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),zbot,pbot,th_bot !,'bot1'
    if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i)
      write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zwat(nzs+1),0.d0,ths(zo(i)) !,'watd1'
    end if
    if(ldeep ) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),zdeep,pdeep,th_deep !,'dp1'
  else if (dcf<=0.d0) then ! Saturated infiltration model
    write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),ztop,ptop,1.d0 !,'top2'
    if(loxlo) then
      if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),ztop-zinc*float(mrpx(1)),p(mrpx(1)),thz(mrpx(1)) !,'max3'
      do m=mfrst, mlast-1 ! Check for wt between min & max 4/25/2012 RLB
        if(wtptr(m)==1 .or. wtptr(m)==2) then
          write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(m),0.d0,ths(zo(i))
        end if
      end do
      if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),ztop-zinc*float(mrpn(1)),p(mrpn(1)),thz(mrpn(1)) !,'min3'
    else
      if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),ztop-zinc*float(mrpn(1)),p(mrpn(1)),thz(mrpn(1)) !,'min4'
      do m=mfrst, mlast-1 ! Check for wt between min & max 4/25/2012 RLB
        if(wtptr(m)==1 .or. wtptr(m)==2) then
          write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(m),0.d0,ths(zo(i))
        end if
      end do
      if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),ztop-zinc*float(mrpx(1)),p(mrpx(1)),thz(mrpx(1)) !,'max4'
    end if
    do m=mlast, nzs ! find lower water table
      if(wtptr(m)==1 .or. wtptr(m)==2) then !
        if (zwat(m) > (zbot*1.00001)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & ix(i),jy(i),zwat(m),0.d0,ths(zo(i)) !,'watu2'
      end if    
    end do
    write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),zbot,pbot,1.d0 !,'bot2'
    if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i) 
      write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zwat(nzs+1),0.d0,ths(zo(i)) !,'watd2'
    end if
    if(ldeep) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
      & ix(i),jy(i),zdeep,pdeep,th_deep !,'dp2'
  end if
end if
return
11 continue
  write(*,*) 'Error writing ijz output file, at step ',jf
  write(ulog,*) 'Error writing ijz output file, at step ',jf
  close (ulog)
  stop '11 in subroutine svijz()'
end subroutine svijz
