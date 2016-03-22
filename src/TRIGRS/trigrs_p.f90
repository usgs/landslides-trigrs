!  Program to compute pore-pressure response and 
!  factor of safety for saturated and unsaturated infiltration.
!  by Rex L. Baum and W.Z. Savage, USGS
!
! MPI calls & MPI-related changes by M. Alvioli, November 2014  
!
program trigrs_mpi
  use input_file_defs; use input_vars
  use grids; use model_vars
  use mpi
  use modules_p
  use partial_arrays_p
  implicit none
  integer, parameter::ulen=25 !Added 12/7/2010 RLB
  integer:: grd
  integer:: i,j,k,imx1,mnd !,m ! "m" removed 1 Feb 2013, RLB
  integer:: nodata,sctr, umax,patlen !added umax, 12/6/2010, RLB
  integer:: ncol,nrow,u(ulen),maxzo,ncc,nccs
  integer:: time_incr_ctr ! added 26 June 2013, RLB 
  real::x1,mnzmx,mndep !,per_dur_min,per_dur
  real::outp_incr_min,outp_incr
  real (double)::newdep,dh ! added 17Nov2014, RLB
  character (len=1):: tb
  character (len=255):: outfil,infil
  character (len=14):: fminfil='TRfs_min_'
  character (len=14):: zfminfil='TRz_at_fs_min_'
  character (len=14):: pminfil='TRp_at_fs_min_'
  character (len=8):: wtabfil='TRwater_'
  character (len=18):: profil='TRlist_z_p_fs_'
  character (len=14):: header(6)
  character (len=13):: ncvfil='TRnon_convrg_'  	
  character (len=8):: date
  character (len=10):: time
  character (len=4):: stp 
  character (len=31):: scratch,irfil
  character (len=7):: vrsn
  character (len=11):: bldate
  character (len=2)::pid(3)
  logical :: lwarn !, lwarn2
  integer myrank,isize,ierr
  Call MPI_INIT(ierr)   ! MPI initialization calls
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  ! first executable statement ............	
  call date_and_time(date,time)
  test=-9999.D0; test1=-9999.
  u=(/11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35/)
  pid=(/'TI','GM','TR'/)
  pi=3.141592653589793
  dg2rad=pi/180.D0
  vrsn='2.0.10t'; bldate='12 Jan 2016'
  smt=0.1d0; lard=12.d0 ! test values for early-time (moved 2/15/2012, smt corrected 3/5/2015, RLB)
  mnd=6 !default value assumed if no integer grid is read.
  fminfil=adjustl(fminfil);zfminfil=adjustl(zfminfil)
  profil=adjustl(profil); pminfil=adjustl(pminfil); 
  if(myrank.eq.0) then
     write (*,*) ''
     write (*,*) 'TRIGRS: Transient Rainfall Infiltration'
     write (*,*) 'and Grid-based Regional Slope-Stability'
     write (*,*) '               Analysis'
     write (*,*) '       Version ', vrsn,', ',bldate
     write (*,*) '  By Rex L. Baum and William Z. Savage'
     write (*,*) '       U.S. Geological Survey'
     write (*,*) '-----------------------------------------'
     write (*,*) 'Portions of this program include material'
     write (*,*) '           copyrighted (C) by'
     write (*,*) '      Absoft Corporation 1988-2010.'
     write (*,*) ''
  endif
  tb=char(9)
  !  open log file ! 14 Feb 2013 Added adjustl() statements for compatability with other compilers
  outfil='TrigrsLog.txt'; outfil=adjustl(outfil)
  if(myrank.eq.0) then
     open (u(19),file=trim(outfil),status='unknown',err=410)
     write (u(19),*) ''
     write (u(19),*) 'Starting TRIGRS ', vrsn,' ',bldate
     write (u(19),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
     write (u(19),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  endif
  !  read initialization file
  call trini_p(u(19),u(7),dg2rad)
  ! determine grid size parameters RLB 4/18/2011
  patlen=scan(elevfil,'/\',.true.) ! find end of folder name
  elfoldr=elevfil(1:patlen) ! path to elevation grid
  ans=.false.
  do i=1,3
     elfoldr=adjustl(elfoldr)
     infil=trim(elfoldr)//pid(i)//'grid_size.txt'
     infil=adjustl(infil)
     if(myrank.eq.0) then
        inquire (file=trim(infil),exist=ans)
        write(*,*) trim(infil), ans
     endif
     if(ans) exit
  end do
  Call MPI_BCAST(ans,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  if(ans) then
     if(myrank.eq.0) then
        open (u(22),file=trim(infil),status='unknown',err=420)
        read (u(22),*) heading
        read (u(22),*) imax,row,col,nwf
        close (u(22))
     endif
     Call MPI_BCAST(imax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     Call MPI_BCAST(row,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     Call MPI_BCAST(col,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     Call MPI_BCAST(nwf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  else
     infil=elevfil; infil=adjustl(infil)
     call ssizgrd_p(row,col,celsiz,nodat,imax,u(12),infil,header,u(19))
     outfil=trim(elfoldr)//pid(3)//'grid_size.txt'
     outfil=adjustl(outfil)
     if(myrank.eq.0) then
        open (u(22),file=trim(outfil),status='unknown',err=410)
        write (u(22),*) 'imax      row      col      nwf'
     endif
     nwf=1 ! dsctr is computed by TopoIndex; dsctr=1 is default value for no runoff routing.
     if(myrank.eq.0) then
        write (u(22),*) imax,row,col,nwf
        write (u(22),*) ''
        close (u(22))
     endif
  end if
  ! 03/2015 set up partial arrays for parallel execution
  call partial_p(imax)
  if(myrank.eq.0) then
     write(u(19),*) 'Grid size parameters from ', trim(infil)
     write (u(19),*) heading
     write (u(19),*) imax,row,col,nwf
  endif
  ! Allocate & initialize arrays needed for runoff routing
  grd=row*col
  imx1=imax
  allocate (pf2(grd),indx(imax),nxt(imax))
  allocate (dsctr(imax+1),slo(imax))
  allocate (pf1(grd),rizero(imax))
  allocate (ri(imax),rik(imax*nper),ro(imax))
  allocate (rikzero(imax),temp(col),itemp(col))
  allocate (depth(imax),zmax(imax))
  allocate (zo(imax),ir(imax),tfg(imax))
  allocate (elev(imax)) ! added 4/21/2010
  elev=0.
  pf2=0
  indx=0
  nxt=0
  dsctr=0
  zo=1
  pf1=0.
  slo=0.
  rizero=0.
  ri=0.
  rik=0.
  ro=0.
  rikzero=0.
  temp=0.;itemp=0
  ir=0.
  depth=0.
  zmax=0.
  ! Choose file extension for grid files, Added 4/14/2010
  grxt='.txt'
  if(lasc) grxt='.asc'
  ! *****************************************************************
  !  read gridded data from GIS
  if(myrank.eq.0) then
     write (*,*) 'Reading input grids'
     write(u(19),*) 'Input file name,            Cell count'
  endif
  !  read slope angles
  call srdgrd_p(grd,col,ncol,nrow,celsiz,nodat,&
       &slo,pf1,sctr,imax,temp,u(1),slofil,param,header,u(19))
  if(myrank.eq.0) then
     write(u(19),*) 'Slope angle grid'
     write(u(19),*) trim(slofil),sctr,' data cells'
     if(sctr/=imax .or. ncol/=col .or. nrow/=row) then
        write(*,*) 'Grid mismatch: ', trim(slofil)
        write(*,*) 'Check slope grid and (or) initialization file against elevation grid.'
        write(u(19),*) 'Grid mismatch: ', trim(slofil)
        write(u(19),*) 'Check slope grid and (or) initialization file against elevation grid.'
     end if
  endif
  slo=slo*dg2rad ! convert slope angles to radians
  !  read property zone numbers, zo
  if(nzon==1) then
     zo=1 ! if only one zone, all values of zone grid equal 1.
     if(myrank.eq.0) then
        write(*,*) 'One property zone, no grid required!'
        write(u(19),*) 'One property zone, no grid required!'
     endif
     parami=param ! added 7/29/2008 RLB
  else
     call irdgrd_p(grd,col,ncol,nrow,celsiz,nodata,mnd,&
          &zo,pf2,sctr,imax,itemp,u(15),zonfil,parami,header,u(19))
     if(myrank.eq.0) then
        write(u(19),*) 'Property zone grid'
        write(u(19),*) trim(zonfil),sctr,' data cells'
     endif
     if(sctr/=imax .or. ncol/=col .or. nrow/=row) then
        if(myrank.eq.0) then
           write (*,*) 'Grid mismatch in TRIGRS main program ',trim(zonfil)
           write (*,*) 'Correct property-zone grid and/or initializtion file.' 
           write (u(19),*) 'Grid mismatch in TRIGRS main program ',trim(zonfil)
           write (u(19),*) 'Correct property-zone grid and/or initializtion file.' 
           close(u(19))
        endif
        call MPI_FINALIZE(ierr)
        stop '-1'
     end if
     maxzo=maxval(zo)
     if (maxzo/=nzon) then
        if(myrank.eq.0) then
           write (*,*) 'Maximum zone number does not equal number of property zones!'
           write (*,*) 'Correct property-zone grid and/or initializtion file.' 
           write (u(19),*) 'Maximum zone number does not equal number of property zones!'
           write (u(19),*) 'Correct property-zone grid and/or initializtion file.' 
           close(u(19))
        endif
        call MPI_FINALIZE(ierr)
        stop '-1'
     end if
  end if
  ! *********************     	
  !  read background infiltration rate, Isub0 
  if (crizero.lt.0) then 
     call srdgrd_p(grd,col,ncol,nrow,celsiz,nodat,&
          &rizero,pf1,sctr,imax,temp,u(16),rizerofil,param,header,u(19))
     if(myrank.eq.0) then
        write(u(19),*) 'Background infiltration rate grid'
        write(u(19),*) trim(rizerofil),sctr,' data cells'
        if(sctr/=imax .or. ncol/=col .or. nrow/=row) write (u(19),*) 'Grid mismatch ',trim(rizerofil)
     endif
  else
     rizero=crizero
  end if
  !  read initial depth to water table, 
  if (dep.lt.0) then  
     call srdgrd_p(grd,col,ncol,nrow,celsiz,nodat,&
          &depth,pf1,sctr,imax,temp,u(10),depfil,param,header,u(19))
     if(myrank.eq.0) then
        write(u(19),*) 'Initial water-table depth grid'
        write(u(19),*) trim(depfil),sctr,' data cells'
        if(sctr/=imax .or. ncol/=col .or. nrow/=row) write (u(19),*) 'Grid mismatch ',trim(depfil)
     endif
  else
     depth=dep
  end if
  !  read depth to base of potential slide, zmax
  if (czmax.lt.0) then 
     call srdgrd_p(grd,col,ncol,nrow,celsiz,nodat,&
          &zmax,pf1,sctr,imax,temp,u(11),zfil,param,header,u(19))
     if(myrank.eq.0) then
        write(u(19),*) 'Maximum depth grid'
        write(u(19),*) trim(zfil),sctr,' data cells'
        if(sctr/=imax .or. ncol/=col .or. nrow/=row) write (u(19),*) 'Grid mismatch ',trim(zfil)
     endif
  else
     zmax=czmax
  end if
  ! Trap error conditions for zmin values.
  mndep=minval(depth) !Added 28 Jan 2013, RLB 
  mnzmx=minval(zmax) !Added 28 Jan 2013, RLB 
  if(zmin>mnzmx .or. zmin>mndep) zmin=0.
  !  read digital elevations, elev !! added 2/24/ 2010, unit number changed to 12 12/6/2010
  call srdgrd_p(grd,col,ncol,nrow,celsiz,nodat,&
       &elev,pf1,sctr,imax,temp,u(12),elevfil,param,header,u(19))
  if(myrank.eq.0) then
     write(u(19),*) 'Elevation grid'
     write(u(19),*) trim(elevfil),sctr,' data cells'
     if(sctr/=imax .or. ncol/=col .or. nrow/=row) write (u(19),*) 'Grid mismatch ',trim(zfil)    
     ! *****************************************************************
     write(u(19),*) '---------------******---------------'
  endif
  ! test and adjust (if necessary) steady background infiltration rates
  call steady_p(sumex,u(19),imx1)
  ! conduct runoff routing and adjust transient infiltration rates
  call rnoff_p(grd,sumex,imx1,celsiz,param,parami,nodat,&
       &nodata,mnd,sctr,ncol,nrow,header,test1,u,ulen)
  ! Deallocate arrays that are no longer needed
  deallocate (ri,ro,wf)
  deallocate (pf2,indx,nxt,dsctr,dsc)
  if(myrank.eq.0) then
     write(u(19),*) '---------------******---------------'
     write(u(19),*) 'Input file name,          Cell count'
  endif
  ! *****************************************************************
  ! compute pore pressure distributions for either fully saturated or
  ! partially saturated conditions.
  ! Partially saturated zone overlies saturated zone
  ! Allocate and initialize new arrays
!  allocate (fsmin(imax*nout),pmin(imax*nout),zfmin(imax*nout))
!  allocate (fsmin_0(imax*nout),pmin_0(imax*nout),zfmin_0(imax*nout))
  allocate (fsmin(isct(myrank)*nout))
  if(myrank.eq.0) allocate(fsmin_0(imax*nout))
  allocate (pmin(isct(myrank)*nout))
  if(myrank.eq.0) allocate(pmin_0(imax*nout))
  allocate (zfmin(isct(myrank)*nout))
  if(myrank.eq.0) allocate(zfmin_0(imax*nout))
  allocate (p(nzs+1),ptran(nzs+1),pzero(nzs+1),bline(nzs+1))
  allocate (fc(nzs+1),fw(nzs+1),thz(nzs+1),kz(nzs+1),trz(nzs+1))
  allocate (nvu(imax),nv(imax),uwsp(nzs+1),gs(nzon))
  allocate (nvu_0(imax),nv_0(imax))
  allocate (chi(nzs+1))
  if(outp(1)) allocate(wtab(imax*nout),wtab_0(imax*nout))
  if(flag<0 .or. outp(1)) then  !Added 17 Nov 2014, RLB
     allocate(p3d(imax*nout,nzs+1),p3d_0(imax*nout,nzs+1))
     allocate(dh3d(imax*nout),newdep3d(imax*nout))  ! moved 14 Jan 2015, RLB
     allocate(dh3d_0(imax*nout),newdep3d_0(imax*nout))
     dh3d=0.d0;newdep3d=0.d0
     dh3d_0=0.d0;newdep3d_0=0.d0
     p3d=0.d0; p3d_0=0.d0
  end if
  if(flag==-1) then
     allocate(fs3d(imax*nout,nzs+1),fs3d_0(imax*nout,nzs+1))
     fs3d=0.d0; fs3d_0=0.d0
  end if
  if(flag==-2) then !Added 17 Nov 2014, RLB
     allocate(pzero3d(imax,nzs+1),ptran3d(imax*nout,nzs+1),fs3d(imax*nout,nzs+1))
     allocate(pzero3d_0(imax,nzs+1),ptran3d_0(imax*nout,nzs+1),fs3d_0(imax*nout,nzs+1))
     pzero3d  =0.d0;ptran3d  =0.d0;fs3d  =0.d0
     pzero3d_0=0.d0;ptran3d_0=0.d0;fs3d_0=0.d0
  end if
  if(flag==-3 .and. unsat0) then !Added 17 Nov 2014, RLB
     allocate(fs3d(imax*nout,nzs+1),th3d(imax*nout,nzs+1))
     allocate(fs3d_0(imax*nout,nzs+1),th3d_0(imax*nout,nzs+1))
     fs3d  =0.d0;th3d  =0.d0
     fs3d_0=0.d0;th3d_0=0.d0
  else if(flag==-3.and. .not. unsat0) then
     allocate(fs3d(imax*nout,nzs+1))
     allocate(fs3d_0(imax*nout,nzs+1))
     fs3d=0.d0; fs3d_0=0.d0
  end if
  if (flag<= -4 .or. outp(1)) then ! moved 5/3/2010, flag =-4,...,-9 produces ijz or xyz output 12/22/2010 RLB
     allocate(ix(imax),jy(imax))
     ix=0;jy=0
     if(flag/=-3) allocate(th3d(imax*nout,nzs+1))
     if(flag/=-3) allocate(th3d_0(imax*nout,nzs+1))
     th3d=0.d0; th3d_0=0.d0
  end if
  fsmin=0; zfmin=0.; pmin=0.
  if(myrank.eq.0) then
     fsmin_0=0.; zfmin_0=0.; pmin_0=0.  
  endif
  if(outp(1)) then
     wtab=0.
     wtab_0=0.
  endif
  p=0.
  ptran=0.
  pzero=0.
  bline=0.
  fc=0.
  fw=0.
  nv=0; nv_0=0
  nvu=0; nvu_0=0
  ! determine number of time steps needed	  
  kper=nper
  if (t>capt(nper+1)) then
     kper=nper+1 
  else
     do k=1,nper ! find the period that contains t
        if(t>=capt(k) .and. t<=capt(k+1)) kper=k
     end do
  end if
  if (tx<1) tx=1
  ! compute time-step duration for unsaturated models, modified for automatic adjustment, 26 June 2013 RLB        
  nts=kper*tx ! number of time-steps from 0 to t
  tns=float(nts)
  tmin=0.
  tmax=t 
  tinc=(tmax-tmin)/tns !timestep duration for unsaturated models
  if(myrank.eq.0) then
     write (u(19),*) 'Initial size of time-steps ', tinc
     write (*,*) 'Initial size of time-steps ', tinc
  endif
  outp_incr_min=tmax! ;per_dur_min=tmax
  if(unsat0)then
     do  k=1,nout-1 ! minimum elapsed time between output 
        if(nout>1) then
           outp_incr=tsav(k+1)-tsav(k)
           if(outp_incr<outp_incr_min) outp_incr_min=outp_incr
        else
           outp_incr=tsav(1)
        end if
     end do
     if(myrank.eq.0) write(*,*)  'outp_incr_min ', outp_incr_min 
     do while (tinc>outp_incr_min)
        tx=tx+1
        nts=kper*tx ! number of time-steps from 0 to t
        tns=float(nts)
        tinc=(tmax-tmin)/tns 
     end do
     if(myrank.eq.0) then
        write (u(19),*) 'Adjusted size of time-steps ', tinc
        write (*,*) 'Adjusted size of time-steps, tx ', tinc, tx
     endif
  endif
  ! compute timestep duration for each period--saturated models, Added 25 June 2013, RLB 
  allocate(tinc_sat(nts+1)); tinc_sat=0.d0 !Revised to eliminate array-boundary errors in  satinf and satfn, 2 Nov 2013, RLB 
  time_incr_ctr=0
  do k=1,kper 
     do i=1,tx
        time_incr_ctr=time_incr_ctr+1
        if(t>=capt(k+1)) then ! Revised 3/6/2015
          tinc_sat(time_incr_ctr)=(capt(k+1)-capt(k))/float(tx)
        else 
          tinc_sat(time_incr_ctr)=(t-capt(k))/float(tx)
        end if
        if(myrank==0 .and. k==kper) write(u(19),*)&
        & 'time_incr_ctr, tinc_sat(time_incr_ctr) ',time_incr_ctr, tinc_sat(time_incr_ctr) ! Revised 3/3/2015, RLB
     end do
  end do
  ! compute output pointers	
  allocate(jsav(nts+1))  
  jsav=0
  if(myrank.eq.0) then
     write (u(19),*) '******** Output times ********'
     write (u(19),*) 'number, timestep #,  time'
     write (*,*) '******** Output times ********'
     write (*,*) 'number, timestep #,  time'
  endif
  lwarn=.false. !; lwarn2=.false.
  do k=1,nout
     ts=tmin
     do j=1,nts
        if(unsat0) then
           if(tsav(k)>=ts .and. tsav(k)<ts+tinc) then
              if(tsav(k)/=ts) lwarn=.true.
              jsav(j)=k
              ksav(k)=j
              tsav(k)=ts
              exit
           else if(tsav(k)>=tmax) then
              jsav(nts+1)=k
              ksav(k)=nts+1
              tsav(k)=tmax
              exit
           end if
           ts=ts+tinc
        else
           if(tsav(k)>=ts .and. tsav(k)<ts+tinc_sat(j)) then
              if(tsav(k)/=ts) lwarn=.true.
              jsav(j)=k
              ksav(k)=j
              tsav(k)=ts
              exit
           else if(tsav(k)>=tmax) then
              jsav(nts+1)=k
              ksav(k)=nts+1
              tsav(k)=tmax
              exit
           end if
           ts=ts+tinc_sat(j)
        endif
     end do
     !            if(unsat0 .and. k>1 . and. tsav(k)==tsav(k-1)) lwarn2 ! Check for matching output times.
     if(lwarn.and.myrank.eq.0) then
        write(u(19),*) 'One or more specified output times unavailable, '
        write(u(19),*) 'Nearest available time substituted.'
        write(*,*) 'One or more specified output times unavailable, '
        write(*,*) 'Nearest available time substituted.'
     end if
     if(myrank.eq.0) then
        write(u(19),*) k,ksav(k),tsav(k)
        write(*,*) k,ksav(k),tsav(k)
     endif
  end do
  ! allocate and initialize additonal model arrays	  
  allocate (r(nmax),q(kper),qtime(2*nts+1),qb(nts+1),tcap(nts+2)) ! qb dimension changed to nts+1 29 Jan 2013, RLB 
  allocate(qts(nts+1)) ! qts assigns q, transient surface flux, to timesteps
  if(outp(7)) then
     allocate(rik1(imax*(nts+1))) ! corrected 6/30/2008, RLB
     rik1=0.
  end if
  eps=1.0e-18
  nmax2=0
  nmin=1+nmax; nmn=1+mmax
  ncc=0;nccs=0
  tis=tiny(x1)
  r=0.; tcap=0; qts=0. ! Added 21 Feb 2013, RLB 
  if(myrank.eq.0) write(*,*) 'Starting computations of pressure head and factor of safety'
  ! prepare file, header, and arrays for generating list files, Added 4/21/2010, revised 5/3/2010, 12/6/2010
  if(flag < 0) then
     umax=maxval(u) ! added 12/6/2010, RLB
     uijz(1)=umax+1
     zmn(1)=minval(elev);zmx(1)=maxval(elev) ! added 9/14/2011, RLB
     if(myrank.eq.0) then
        write(u(19),*) 'DEM minimum & maximum elevations:', zmn(1), zmx(1) ! added 2/14/2012 RLB
        call prpijz(u(2),u(19),profil,col,row,header,vrsn) !Revised 12/23/2010
     endif
  end if
  ! ---------------------------------------------------------------
  if(myrank.eq.0) then
     call date_and_time(date,time)
     write (u(19),*) 'Starting TRIGRS ', vrsn,' ',bldate
     write (u(19),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
     write (u(19),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  endif
  ! ---------------------------------------------------------------
  if(unsat0) then
     do j=1,nzon ! compute specific gravity of solids (gs(j)) from saturated unit weight of soil (uws(j))
        gs(j)=((uws(j)/uww)-ths(j))/(1-ths(j))
     end do
     if(mmax.lt.0) then ! infinite depth model
        mmax=20; nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
        if(myrank.eq.0) then
           write(*,*) 'Calling unsaturated infinite-depth model'
           !  4/14/2010 RLB added logging of main subroutines handling infiltration 	    
           write(u(19),*) 'Unsaturated infinite-depth model, unsinf()'
        endif
        call unsinf_p(imx1,u(19),u(2),ncc,nccs)
     else ! finite depth model
        if(myrank.eq.0) then
           write(*,*) 'Calling unsaturated finite-depth model'
           write(u(19),*) 'Unsaturated finite-depth model, unsfin()'
        endif
        call unsfin_p(imx1,u(19),u(2),ncc,nccs)
     end if
  else ! Saturated zone extends to ground surface
     if(myrank.eq.0) write(*,*) 'Ignoring unsaturated zone'
     if(tx==1 .and. nout==1) then
        deallocate (rizero,ir)
        !  compute pore-pressure distributions and factor of safety
        outfil=trim(folder)//trim(profil)//trim(suffix)//'.txt'
        nccs=0
        if(mmax.lt.0) then ! infinite depth model
           mmax=20; nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
           if(myrank.eq.0) then
              write(*,*) 'Calling saturated infinite-depth model'
              write(u(19),*) 'Saturated infinite-depth model, iverson()'              
              write(*,*) " - - - - - - WARNING: parallel version unavailable - - - - -"
              write(u(19),*) " - - - - - - WARNING: parallel version unavailable - - - - -"
           endif
           call iverson(imx1,u(2),outfil,u(19))
        else ! finite depth model
           if(myrank.eq.0) then
              write(*,*) 'Calling saturated finite-depth model'
              write(u(19),*) 'Saturated finite-depth model, savage()'
              write(*,*) " - - - - - - WARNING: parallel version unavailable - - - - -"
              write(u(19),*) " - - - - - - WARNING: parallel version unavailable - - - - -"
           endif
           nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
           call savage(imx1,u(2),outfil,u(19),nccs)
        end if
     else
        if(mmax.gt.0) then	
           if(myrank.eq.0) then
              write(*,*) 'Calling multistep saturated finite-depth model'
              write(u(19),*) 'Multistep saturated finite-depth model, satfin()'
           endif
           nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
           call satfin_p(imx1,u(19),u(2),nccs)
        else
           if(myrank.eq.0) then
              write(*,*) 'Calling multistep saturated infinite-depth model'
              write(u(19),*) 'Multistep saturated infinite-depth model, satinf()'
           endif
           mmax=20; nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
           call satinf_p(imx1,u(19),u(2),nccs)
        end if
     end if
  end if
  ! *****************************************************************
  !  write output grid files 
  ! 4/14/2010 added option to let file extension be either ".asc" or ".txt"
101 continue
  if(myrank.ne.0) goto 1999
  call date_and_time(date,time)
  write (u(19),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
  write (u(19),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  write(*,*) 'Saving results'
  ti=tiny(param(1)) ! Changed from param(m) to param(1), 28 Jan 2013, RLB
  do j=1,nout
     write(stp,'(i4)') j
     stp=adjustl(stp)
     if (outp(3)) then ! minimum factor of safety
        tfg=0.
        do i=1,imx1
           tfg(i)=fsmin_0(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(fminfil)//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,u(4),test1,param,u(19),outfil,ti,header)
     end if
     if (outp(4)) then ! depth of minimum factor of safety
        tfg=0.
        do i=1,imx1
           tfg(i)=zfmin_0(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(zfminfil)//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,u(5),test1,param,u(19),outfil,ti,header)
     end if
     if (outp(5)) then ! pressure head at depth of minimum factor of safety
        tfg=0.
        do i=1,imx1
           tfg(i)=pmin_0(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(pminfil)//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,u(6),test1,param,u(19),outfil,ti,header)
     end if
     if (flag<=-4 .or. outp(1)) then ! Moved from computational subroutines 17Nov2014, RLB
        if (flag>=-6) then
           do i=1,imx1
              newdep=newdep3d(i+(j-1)*imax)
              dh=dh3d(i+(j-1)*imax)
              do k=1,nzs+1
                 ! add code to handle saturated cells ! dcf=0.d0; chi=1.d0 
                 p(k)=p3d(i+(j-1)*imax,k)
                 thz(k)=th3d(i+(j-1)*imax,k)
              end do
              call svijz(i,j,dh,newdep,u(19))
           end do
        endif
        if (flag<=-7 .and. flag>=-9) then   ! Added 2/10/2012
           do i=1,imx1
              newdep=newdep3d(i+(j-1)*imax)
              dh=dh3d(i+(j-1)*imax)
              do k=1,nzs+1
                 p(k)=p3d(i+(j-1)*imax,k)
                 thz(k)=th3d(i+(j-1)*imax,k)
              end do
              call svxmdv(i,j,dh,newdep,u(19))
           end do
        endif
     end if
     if (outp(1)) then ! computed water table
        tfg=0.
        do i=1,imx1
           tfg(i)=wtab(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(wtabfil)//el_or_dep//'_'//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,u(6),test1,param,u(19),outfil,ti,header)
     end if
  end do
  if(flag<=-1 .and. flag>=-3) then
     !	      write(*,*)'j,tsav(j)',j,tsav(j)
     call svlist(u(2),tsav(j),j)
  end if  
  if (outp(7) .and. unsat0) then ! incremental basal flux, unsaturated zone
     do j=1,nts
        ir=0.
        do i=1,imx1
           ir(i)=rik1(i+(j-1)*imax)
        end do
        ti=tiny(param(1))
        write(scratch,'(i6)') j
        scratch=adjustl(scratch)
        irfil='TRunszfluxTS'//trim(scratch)//trim(suffix)//grxt
        irfil=adjustl(irfil)
        outfil=trim(folder)//trim(irfil)
        call ssvgrd(ir,imax,pf1,row,col,u(9),test1,param,u(19),outfil,ti,header)
     end do
  end if
  if (ncc>0) then ! non-convergent cells, unsaturated zone (12/6/2010, changed unit # from 7 to 14)
     outfil=trim(folder)//ncvfil//'UZ_'//trim(suffix)//grxt
     call isvgrd(nvu,imax,pf1,row,col,u(14),test,test,mnd,parami,u(19),outfil,ti,header)
  end if
  if (nccs>0) then ! non-convergent cells, saturated zone
     outfil=trim(folder)//ncvfil//'SZ_'//trim(suffix)//grxt
     call isvgrd(nv,imax,pf1,row,col,u(14),test,test,mnd,parami,u(19),outfil,ti,header)
  end if
  write (*,*) 'TRIGRS finished!'
  write (u(19),*) 'TRIGRS finished normally'
  if(flag<=-1) close(u(2)) ! moved from subroutines 1 Dec 2011 RLB
  call date_and_time(date,time)
  write (u(19),*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4)
  write (u(19),*) 'Time ',time(1:2),':',time(3:4),':',time(5:6)
  close (u(19))
  goto 1999
  stop '0'
  ! Error reporting 	
410 continue
  write (*,*) 'Error opening output file in TRIGRS main program'
  write (*,*) '--> ',outfil
  write (*,*) 'Check file path and status'
  write (u(19),*) 'Error opening output file in TRIGRS main program'
  write (u(19),*) '--> ',outfil
  write (u(19),*) 'Check file path and status'
  call MPI_FINALIZE(ierr)
  stop '410'
420 continue
  write (*,*) 'Error opening input file in TRIGRS main program'
  write (*,*) '--> ',infil
  write (*,*) 'Check file path and status'
  write (u(19),*) 'Error opening input file in TRIGRS main program'
  write (u(19),*) '--> ',infil
  write (u(19),*) 'Check file path and status'
  call MPI_FINALIZE(ierr)
  stop '420'
  !
  1999 continue
  call MPI_FINALIZE(ierr)
  stop
2000 continue
  if(myrank.eq.0) then
     print *
     print *," - - - - - - - - - parallel version unavailable - - - - - - - -"
     print *
  endif
  call MPI_FINALIZE(ierr)
  stop
end program
