!
! MPI calls & MPI-related changes by M. Alvioli, November 2014  
!
subroutine trini_p(ulog,uini,dg2rad)
  ! reads initialization file, R.L. Baum, USGS, latest revision 29 Mar 2013 
  use input_file_defs; use input_vars; use model_vars, only: smt
  use mpi
  implicit none
  integer:: i,j,iz,linct
  integer, intent(in):: ulog,uini
  real, intent(in):: dg2rad
  real :: tstar,tdif
  logical :: ltdif
  character (len=31):: scratch
  integer myrank,isize,ierr
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  init='tr_in.txt'
  init=adjustl(init)
  if(myrank.eq.0) then
     inquire (file=trim(init),exist=ans)
     if(ans) then
        open (uini,file=trim(init),status='old',err=201)
        write (*,*) 'Opening default initialization file'
     else
        write (*,*) 'Cannot locate default initialization file, <tr_in.txt>'
        write (*,*) 'Type name of initialization file and'
        write (*,*) 'press RETURN to continue'
        read (*,'(a)') init
        init=adjustl(init)
        open (uini,file=trim(init),status='old',err=201)
     end if
     write (ulog,*) 'initialization file -->',init
     write (ulog,*) '-- LISTING OF INITIALIZATION FILE --'	
     ! write copy of data to log file
     linct=1	
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) title; linct=linct+1
     title=adjustl(title)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(title)
     write (*,*) title
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) tx,nmax,mmax,nzon; linct=linct+1
     write (ulog,*) trim(heading)
     if(nmax<2) nmax=2 ! set minimum value for nmax
     write (ulog,*) tx,nmax,mmax,nzon
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) nzs,zmin,uww,nper,t
     linct=linct+1
     write (ulog,*) trim(heading)
     write (ulog,*) nzs,zmin,uww,nper,t
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) czmax,dep,crizero,slomin,slomax; linct=linct+1
     write (ulog,*) trim(heading)
     write (ulog,*) czmax,dep,crizero,slomin,slomax
  endif
  Call MPI_BCAST(tx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(nmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(mmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(nzon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(nzs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(zmin,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(uww,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(nper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(t,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(czmax,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(dep,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(crizero,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(slomin,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(slomax,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  if(slomax<0. .or. slomax>90.) slomax=90.
  if(slomin<0. .or. slomin>=slomax) slomin=0.
  slomin=slomin*dg2rad ! convert minimum slope angle to radians
  slomax=slomax*dg2rad ! convert maximum slope angle to radians
  ! allocate & read arrays for zone properties and initial conditions  	
  allocate (c(nzon),phi(nzon),uws(nzon),dif(nzon),ks(nzon),ths(nzon),thr(nzon),alp(nzon),unsat(nzon))
  allocate (igcap(nzon)) ! Added 2/15/2012 RLB
  c=0;phi=0
  uws=0;dif=0;ks=0;ths=0;thr=0;alp=0
  unsat=.true.;unsat0=.false.; igcap=.false.; igcapf=.true.
  if(myrank.eq.0) then
     do i=1,nzon
        read (uini,*,err=420) scratch,iz ! property zone number
        linct=linct+1
        scratch=adjustl(scratch)
        read (uini,'(a)',err=420) heading; linct=linct+1
        heading=adjustl(heading) ! zone parameters
        read (uini,*,err=420) c(iz),phi(iz),uws(iz),dif(iz),ks(iz),ths(iz),thr(iz),alp(iz)
        linct=linct+1
        write (ulog,*) trim(scratch),': ',iz
        write (ulog,*) trim(heading)
        write (ulog,*) c(iz),phi(iz),uws(iz),dif(iz),ks(iz),ths(iz),thr(iz),alp(iz)
        if (c(iz)<0. .or. phi(iz)<0. .or. uws(iz)<0. .or.&
             & dif(iz)<0. .or. ks(iz)<0. .or. ths(iz) <0. .or. ths(iz) <0.) then ! check neg. vaules RLB 12/7/2010
           write(*,*) 'Error, negative property value in line ',linct
           write(*,*) 'Edit tr_in.txt and restart program TRIGRS'
           write(ulog,*) 'Error, negative property value in line ',linct
           write(ulog,*) 'Edit tr_in.txt and restart program TRIGRS'
           stop 'Negative property in trini()'
        end if
        if(ths(iz)>0.95) then ! set default value for saturated water content, Added 22 Aug 2013, RLB  
           write(*,*)'Error, Theta-sat. > 0.95 in property zone',iz
           write(*,*)'Theta-sat. value of 0.45 will be used for cells in zone', iz,'.'
           write(ulog,*)'Error, Theta-sat. > 0.95 in property zone',iz
           write(ulog,*)'Theta-sat. value of 0.45 will be used for cells in zone', iz,'.'
        end if
        if(ths(iz)<thr(iz)) then 
           write(*,*)'Error, Theta-resid. > Theta-sat. for property zone',iz
           write(*,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
           write(ulog,*)'Error, Theta-resid. > Theta-sat. for property zone',iz
           write(ulog,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
           unsat(iz)=.false.
        end if
        if(alp(iz)<=0) then
           unsat(iz)=.false.
           write(*,*)'Negative or zero value of Alpha for property zone',iz
           write(*,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
           write(ulog,*)'Negative or zero value of Alpha for property zone',iz
           write(ulog,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
        end if
        if(unsat(iz)) then
           unsat0=.true. ! tracks whether any property zones are unsat.
           write(*,*)'Unsaturated infiltration model selected for cells in zone', iz,'.'
           write(ulog,*)'Unsaturated infiltration model selected for cells in zone', iz,'.'
        end if
     end do
  endif
  Call MPI_BCAST(c,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(phi,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(uws,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(dif,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(ks,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(ths,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(thr,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(alp,nzon,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(unsat,nzon,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(unsat0,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  phi=phi*dg2rad ! convert phi angles to radians
  ! Allocate & read arrays for storm period data	
  allocate (cri(nper),capt(nper+2),rifil(nper)) 
  if(myrank.eq.0) then
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) (cri(j), j=1,nper) ! List of rainfall rates
     linct=linct+1
     write (ulog,*) trim(heading)
     write (ulog,*) (cri(j), j=1,nper)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) (capt(j), j=1,nper+1) ! List of times corresponding to change of rate
     capt(nper+2)=t ! for cases where t>capt(nper+1)
     linct=linct+1
     write (ulog,*) trim(heading)
     write (ulog,*) (capt(j), j=1,nper+1)
     ltdif=.false. ! Test time-step order and error message added 29 Jan 2013, RLB 
     do j=1,nper
        tdif=capt(j+1)-capt(j)
        if(tdif<0.) ltdif=.true. 
     end do
     if(ltdif) goto 424
     !  path names of input files
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) slofil ! File name of slope angle grid (slofil)
     linct=linct+1
     slofil=adjustl(slofil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(slofil)
     ! added 4/21/2010 ******************************************
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) elevfil ! File name of digital elevation grid (elevfil)
     linct=linct+1
     elevfil=adjustl(elevfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(elevfil)
     ! **********************************************************
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) zonfil ! File name of property zone grid (zonfil)
     linct=linct+1
     zonfil=adjustl(zonfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(zonfil)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) zfil ! File name of depth grid (zfil)
     linct=linct+1
     zfil=adjustl(zfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(zfil)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) depfil ! File name of initial depth of water table grid (depfil)
     linct=linct+1
     depfil=adjustl(depfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(depfil)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=420) rizerofil ! File name of initial infiltration rate grid (rizerofil)
     linct=linct+1
     rizerofil=adjustl(rizerofil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(rizerofil)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     write (ulog,*) trim(heading)
     do j=1,nper
        read (uini,'(a)',err=421) rifil(j) ! List of file names of rainfall intensity for each period, (rifil())
        linct=linct+1	
        rifil(j)=adjustl(rifil(j))
        write (ulog,*) trim(rifil(j))
     end do
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=421) nxtfil ! File name of grid of D8 runoff receptor cell numbers (nxtfil)
     linct=linct+1	
     nxtfil=adjustl(nxtfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(nxtfil)
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=421) ndxfil ! File name of list of defining runoff computation order (ndxfil)
     linct=linct+1	
     ndxfil=adjustl(ndxfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(ndxfil)
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=421) dscfil ! File name of list of all runoff receptor cells  (dscfil)
     linct=linct+1	
     dscfil=adjustl(dscfil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(dscfil)
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=421) wffil ! File name of list of runoff weighting factors  (wffil)
     linct=linct+1	
     wffil=adjustl(wffil)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(wffil)
     !  location of output files	
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=421) folder ! Folder where output grid files will be stored  (folder)
     linct=linct+1	
     folder=adjustl(folder)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(folder)
     !  output-file ID code
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,'(a)',err=421) suffix ! Identification code to be added to names of output files (suffix)
     linct=linct+1	
     suffix=adjustl(suffix)
     write (ulog,*) trim(heading)
     write (ulog,*) trim(suffix)
     !  output file selections	
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=421) rodoc ! Save grid files of runoff?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) rodoc
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=421) outp(3) ! Save grid of minimum factor of safety?
     write (ulog,*) trim(heading)
     write (ulog,*) outp(3)
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=421) outp(4) ! Save grid of depth of minimum factor of safety?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) outp(4)
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=421) outp(5) ! Save grid of pore pressure at depth of minimum factor of safety?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) outp(5)
     read (uini,'(a)',err=421) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=422) outp(1), el_or_dep ! Save grid of computed water table depth or elevation for specified time steps?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) outp(1), el_or_dep
     read (uini,'(a)',err=422) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=422) outp(6) ! Save grid files of actual infiltration rate?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) outp(6)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) outp(7) ! Save grid files of unsaturated zone basal flux?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) outp(7)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=423) flag,spcg! Save listing of pressure head and factor of safety ("flag")
     !                              !  & increment (spcg)? ! Added spcg 2/14/2012 RLB
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) flag,spcg! Added spcg 2/14/2012 RLB
     read (uini,'(a)',err=423) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=423) nout ! Number of times to save output grids
     linct=linct+1	
     if(nout<1) nout=1 ! must save at least one time; negative number not allowed
     ! Add code to limit nout to less than some specific number of values?
     write (ulog,*) trim(heading)
     write (ulog,*) nout
  endif
  Call MPI_BCAST(cri,nper,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(capt,nper+2,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(slofil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(elevfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(zonfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(zfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(depfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(rizerofil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  do i=1,nper
     Call MPI_BCAST(rifil(i),255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  enddo
  Call MPI_BCAST(nxtfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(ndxfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(dscfil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(wffil,255,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(folder,224,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(suffix,8,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(rodoc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(outp,8,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(el_or_dep,5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(spcg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(nout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate (tsav(nout),ksav(nout),uijz(nout)) !Added uijz 12/07/2010 RLB
  if(myrank.eq.0) then
     tsav=0.;ksav=0
     read (uini,'(a)',err=420) heading; linct=linct+1 
     heading=adjustl(heading)
     read (uini,*,err=420) (tsav(j), j=1,nout) ! Times of output grids
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) (tsav(j), j=1,nout)
     ! user options  	
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) lskip ! Skip other timesteps?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) lskip
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) lany ! Use analytic solution for fillable porosity?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) lany
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) llus ! Estimate positive pressure head in rising water table zone
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) llus
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) lps0 ! Use psi0=-1/alpha?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) lps0
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) outp(8) ! Log mass balance results?
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) outp(8)
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) flowdir ! Specify flow direction
     linct=linct+1	! 29 Mar 2013, RLB, Added afer each read statement to
     ! correct error in counting line numbers from this point to end of input.
     flowdir=adjustl(flowdir)
     write (ulog,*) trim(heading)
     write (ulog,*) flowdir
     ! added 19 Aug 2009 RLB-------------------
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) bkgrof ! Specify background flux offset 
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) bkgrof
     ! added 14-15 Apr 2010 RLB-------------------
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) lasc ! Specify grid file extension (.asc if true, default is .txt)  
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) lasc
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) lpge0 ! Ignore negative pore pressures in FS (p=0 if true)  
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) lpge0
     ! added 2-15/2012 RLB ----
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) igcapf ! Ignore height of capilary fringe greater than depth 
     !                            ! to water table in selecting infiltration model? 
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) igcapf
     !Parameters for deep pore-pressure estimate in SCOOPS ijz output (Added 20 Apr 2011, RLB): 
     read (uini,'(a)',err=420) heading; linct=linct+1
     heading=adjustl(heading)
     read (uini,*,err=420) deepz,deepwat ! Deep point (floating point) depth below ground surface,  pressure option ('zero' or 'flow')
     linct=linct+1	
     write (ulog,*) trim(heading)
     write (ulog,*) deepz,deepwat
     !-----------------------------------------  	
     close(uini)     
     write (ulog,*) '-- END OF INITIALIZATION DATA --'	
     write (ulog,*) ''
     write (ulog,*) trim(title)
     write (ulog,*) ''
  endif
  Call MPI_BCAST(tsav,nout,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(lskip,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(lany,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(llus,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(lps0,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(flowdir,5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(bkgrof,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(lasc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(lpge0,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(igcapf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  Call MPI_BCAST(deepz,1,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(deepwat,4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  do iz=1,nzon     
     if(alp(iz)>=0 .and. unsat(iz)) then ! Revised 2/15/2012 RLB
        if(igcapf) igcap(iz)=.true. ! Use unsaturated model even for large capillary fringe
        ! Use saturated model for situations where much of interval is small time
        tstar=t*ks(iz)*alp(iz)/(ths(iz)-thr(iz))
        if(tstar < 4.d0*smt) igcap(iz)=.false.
     end if
     if(myrank.eq.0) then
        if (igcap(iz)) then
           write(*,*) '******** Zone ',iz, ' *********'
           write(*,*)'Using unsaturated infiltration model.'
           write(ulog,*) '******** Zone ',iz, ' *********'
           write(ulog,*)'Using unsaturated infiltration model.'
        else if (igcapf) then ! ................................................................ RLB 2/21/2012
           write(*,*) '******** Zone ',iz, ' *********'
           write(*,*)'Using saturated infiltration model to avoid'
           write(*,*)'early-time errors in unsaturated infiltration model.'
           write(ulog,*) '******** Zone ',iz, ' *********'
           write(ulog,*)'Using saturated infiltration model to avoid'
           write(ulog,*)'early-time errors in unsaturated infiltration model.'
        else if (alp(iz)<0) then ! ................................................................ RLB 1/30/2013 
           write(*,*) '******** Zone ',iz, ' *********'
           write(*,*)'Using saturated infiltration model; Alpha<0.'
           write(ulog,*) '******** Zone ',iz, ' *********'
           write(ulog,*)'Using saturated infiltration model; Alpha<0.'
        else
           ! ................................................................ RLB 2/15/2012          
           write(*,*) '******** Zone ',iz, ' *********'
           write(*,*)'Using unsaturated infiltration model.'
           write(*,*) 'Cells where water table is shallower than '
           write(*,*) '           ', 1./alp(iz)
           write(*,*) 'treated as tension saturated--Saturated infiltration model used.'
           write(ulog,*) '******** Zone ',iz, ' *********'
           write(ulog,*)'Using unsaturated infiltration model.'
           write(ulog,*) 'Cells where water table is shallower than '
           write(ulog,*) '           ', 1./alp(iz)
           write(ulog,*) 'treated as tension saturated--Saturated infiltration model used.'
        end if
     endif
  end do
  if(myrank.eq.0) then
     write(*,*) '********  ********  ********  *********'
     write(ulog,*) '********  ********  ********  *********'
  endif
  return
201 continue
  write (*,*) '*** Error opening intialization file in subroutine trini ***'
  write (*,*) '--> ',trim(init)
  write (*,*) 'Check file location and name'
  write (ulog,*) '*** Error opening intialization file in subroutine trini ***'
  write (ulog,*) '--> ',trim(init)
  write (ulog,*) 'Check file location and name'
  write(*,*) 'Press RETURN to exit'
  read*
  stop '201 in trini()'
420 continue
  write (*,*) 'Error reading initialization file in subroutine trini'
  write (*,*) '--> ',trim(init), ' at line ',linct
  write (*,*) 'Check file contents and organization'
  write (ulog,*) 'Error reading initialization file in subroutine trini'
  write (ulog,*) '--> ',trim(init), ' at line ',linct
  write (ulog,*) 'Check file contents and organization'
  call MPI_FINALIZE(ierr)
  stop '420 in trini()'
421 continue
  write (*,*) 'Error reading initialization file in subroutine trini'
  write (*,*) '--> ',trim(init), ' at line ',linct
  write (*,*) 'Check file contents and organization'
  write (*,*) 'Number of file names/place holders for rainfall data'
  write (*,*) 'must equal nper.  List each on a separate line.'
  write (ulog,*) 'Error reading initialization file in subroutine trini'
  write (ulog,*) '--> ',trim(init), ' at line ',linct
  write (ulog,*) 'Check file contents and organization'
  write (ulog,*) 'Number of file names/place holders for rainfall data'
  write (ulog,*) 'must equal nper.  List each on a separate line.'
  call MPI_FINALIZE(ierr)
  stop '421 in trini()'
422 continue
  write (*,*) 'Error reading initialization file in subroutine trini'
  write (*,*) '--> ',trim(init), ' at line ',linct
  write (*,*) 'Check file contents and organization'
  write (*,*) 'Be sure to specify water table elevation or depth'
  write (*,*) 'along with choice or whether or not to save to file.'
  write (ulog,*) 'Error reading initialization file in subroutine trini'
  write (ulog,*) '--> ',trim(init), ' at line ',linct
  write (ulog,*) 'Check file contents and organization'
  write (ulog,*) 'Be sure to specify water table elevation or depth'
  write (ulog,*) 'along with choice or whether or not to save to file.'
  call MPI_FINALIZE(ierr)
  stop '422 in trini()'
423 continue
  write (*,*) 'Error reading initialization file in subroutine trini'
  write (*,*) '--> ',trim(init), ' at line ',linct
  write (*,*) 'Check file contents and organization'
  write (*,*) 'Be sure to specify vertical spacing increment for list file'
  write (*,*) 'on same line as output flag.'
  write (ulog,*) 'Error reading initialization file in subroutine trini'
  write (ulog,*) '--> ',trim(init), ' at line ',linct
  write (ulog,*) 'Check file contents and organization'
  write (ulog,*) 'Be sure to specify vertical spacing increment for list file'
  write (ulog,*) 'on same line as output flag.'
  call MPI_FINALIZE(ierr)
  stop '423 in trini()'
424 continue ! Error message added 29 Jan 2013, RLB 
  write (*,*) 'Subroutine trini reports time steps out of order in file'
  write (*,*) '--> ',trim(init), ' at line ',linct
  write (*,*) 'List time steps in increasing order'
  write (ulog,*) 'Subroutine trini reports time steps out of order in file'
  write (ulog,*) '--> ',trim(init), ' at line ',linct
  write (ulog,*) 'List time steps in increasing order'
  call MPI_FINALIZE(ierr)
  stop '424 in trini()'
  !
end subroutine trini_p

