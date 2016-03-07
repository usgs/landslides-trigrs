	subroutine rnoff(grd,sumex,imx1,celsiz,param,parami,nodat,&
        &nodata,mnd,sctr,ncol,nrow,header,test1,u,unum)
! Runoff routing for TRIGRS        
! By Rex L. Baum, 1 April 2004, latest revision 1 Jul 2011
	use input_file_defs; use input_vars
	use grids
	implicit none
	integer,parameter:: double=kind(1d0)
  	integer:: nodata,sctr,id,next 
	integer:: i,j,l,imx1,roflg,mnd !,m "m" removed 1 Feb 2013, RLB 
	integer:: grd,unum
	integer:: ncol,nrow,u(unum)
	logical:: logi(5)
	character (len=255):: infil,outfil
	character (len=31):: rofil, scratch,irfil
	character (len=14):: header(6)
	real:: rnof,inflx,test1
	real (double):: nodat,sumin,sumro,sumrf,sumex
	real (double):: celsiz,param(6),parami(6),ti
! verify that the runoff routing input files are named in the initialization file	
	logi(5)=.true.
	if(trim(nxtfil)=='no_input') logi(5)=.false.
	if(trim(ndxfil)=='no_input') logi(5)=.false.
	if(trim(dscfil)=='no_input') logi(5)=.false.
	if(trim(wffil)=='no_input') logi(5)=.false.
	ans=logi(5)
	if(logi(5)) then
! .... and that the named files acually exist	
	  inquire (file=trim(nxtfil),exist=logi(1))
	  inquire (file=trim(ndxfil),exist=logi(2))
	  inquire (file=trim(dscfil),exist=logi(3))
	  inquire (file=trim(wffil),exist=logi(4))
	  infil=trim(elfoldr)//'TIgrid_size.txt'
	  infil=adjustl(infil)
	  inquire (file=trim(infil),exist=logi(5)) ! added 2/24/2011, Revised to add folder 21 Feb 2013, RLB
	  do i=1,5
	    if(logi(i)) then
	      continue
	    else
	      ans=.false.
	    end if
	  end do
	end if
	if(ans) then
     	  write(*,*) 'Starting runoff-routing computations'
! Added 2/2/2011, revised 4/1/2011 RLB
  	  if(nwf<imax) then 
  	    nwf=imax*2
  	    write(u(19),*) 'Error, nwf <imax, reset nwf to ',nwf *2
  	    write(u(19),*) 'If errors occur in runoff routing, check "TI_grid_size.txt" for correct value of nwf!'
  	    write(*,*) 'Error, nwf <imax, reset nwf to ',nwf*2 
  	    write(*,*) 'If errors occur in runoff routing, check "TI_grid_size.txt" for correct value of nwf!'
  	  end if
	 allocate (dsc(nwf),wf(nwf))
	 wf=0.; dsc=0
!  read numbers of subjacent cells
	 call irdgrd(grd,col,ncol,nrow,celsiz,nodata,mnd,&
         &nxt,pf2,sctr,imax,itemp,u(17),nxtfil,parami,header,u(19))
    	 write(u(19),*) 'Downslope cell grid'
     	 write(u(19),*) trim(nxtfil),sctr,' data cells'
     	 if(sctr/=imx1) write (u(19),*) 'Grid mismatch ',trim(nxtfil)
! read cell number index to determine order of computation for 
! runoff routing
	 infil=ndxfil
	 open (u(18),file=infil,status='old',err=400)
	 do i=1,imax
 	   read (u(18),*,end=430) j,indx(j)
  	 end do
  	 close (u(18))
! read cell numbers and weighting factors for runoff routing	
	 infil=dscfil
	 open (u(20),file=infil,status='old',err=400)
	 call irdswm(nwf,imax,u(20),nodata,dsc,dsctr,u(19))
	 close (u(20))
	 infil=wffil
	 open (u(21),file=infil,status='old',err=400)
	 call srdswm(nwf,imax,u(21),test1,wf,dsctr,u(19))
	 close (u(21))
!  compute normalized infiltration intensities, Itransient/Ks, rik,
     	 write(u(19),*) ''
 	 write(u(19),*) '*** Runoff-routing computations ***'
	 do j=1,nper
	   sumin=0.d0
	   sumrf=0.d0
	   sumro=0.d0
	   roflg=0
!  read precipitation intensity, I
	   if (cri(j).lt.0) then
	     call srdgrd(grd,col,ncol,nrow,celsiz,nodat,ri,&
     	     &pf1,sctr,imax,temp,u(13),rifil(j),param,header,u(19))
    	   write(u(19),*) 'Precipitation intensity grid ',j
     	   write(u(19),*) trim(rifil(j)),sctr,' data cells'
     	   if(sctr/=imx1) write (u(19),*) 'Grid mismatch ',trim(rifil(j))
  	  else
	    do  i=1,imx1
	      ri(i)=cri(j)
  	    end do
  	  end if
! initialize ro and ir at beginning of each period	
	  do i=1,imx1
	    ir(i)=ks(zo(i))
	    ro(i)=0.
	    sumrf=sumrf+ri(i)
  	  end do
! compute infiltration and runoff at each cell for each time period  
	  do i=1,imx1
	    id=indx(i)
	    next=nxt(id)
	    inflx=ro(id)+ri(id)
! case 1. exfiltration at cells where water table is initially at
! at the ground surface; does not track outflow from cells where
! the water table was initially below the surface and later filled up.
  	    if(depth(id)==0.0 .and. rizero(id)<0.0) then
	      ir(id)=0
	      rik(id+(j-1)*imax)=0.
	      rnof=inflx-rizero(id)
	      ro(id)=rnof
	      do l=dsctr(id), dsctr(id+1)-1
	       if (dsc(l).eq.id) then
	         ro(dsc(l))=ro(dsc(l))+rnof*(wf(l)-1.)   
	       else
	         if (dsc(l)<1) write (*,*) dsc(l)
	         if (dsc(l)>imx1) write (*,*) dsc(l)
	         ro(dsc(l))=ro(dsc(l))+rnof*wf(l)
	       end if
    	      end do	    
	      sumin=sumin+0
! case 2. infiltration, but available water > Ks	    
	    else if (ks(zo(id)).lt.inflx) then
	      rik(id+(j-1)*imax)=1.d0
	      rnof=inflx-ks(zo(id))
	      ro(id)=rnof
	      do l=dsctr(id), dsctr(id+1)-1
	       if (dsc(l).eq.id) then
	         ro(dsc(l))=ro(dsc(l))+rnof*(wf(l)-1.)   
	       else
	         if (dsc(l)<1) write (*,*) dsc(l)
	         if (dsc(l)>imx1) write (*,*) dsc(l)
	         ro(dsc(l))=ro(dsc(l))+rnof*wf(l)
	       end if
    	      end do	    
	      sumin=sumin+ks(zo(id))
	    else
! case 3. available water < Ks
	      ir(id)=inflx
  	      rik(id+(j-1)*imax)=inflx/ks(zo(id))
	      rnof=0.
	      ro(id)=rnof
	      ro(next)=ro(next)+rnof
	      sumin=sumin+inflx
	    end if
  	  end do
! compute total runoff  
  	  do i=1,imx1
	    if (ro(i)>0.) roflg=1
	    if(i==nxt(i)) sumro=sumro+ro(i)
  	  end do
  	  write(u(19),*) ''
  	  write (u(19),*) 'Mass Balance Totals for period ',j
  	  write (u(19),*) 'Precipitation + Exfiltration = &
     &Infiltration + Runoff'
	  write (u(19),*) sumrf+sumex,' : ',sumin+sumro
  	  write (u(19),*) 'Infiltration   Runoff'
	  write (u(19),*) sumin,sumro
  	  write (u(19),*) 'Precipitation   Exfiltration'
	  write (u(19),*) sumrf,sumex
  	  write(u(19),*) '------------------------------------------'
! output ro() arrays for grids that have non-zero runoff
	  if(roflg==1 .and. rodoc) then
! write only grids that are at or before the current time.
	   if(capt(j)<=t) then
	    ti=tiny(param(1))  ! Changed from param(m) to param(1), 28 Jan 2013, RLB
	    write(scratch,'(i6)') j
	    scratch=adjustl(scratch)
	    rofil='TRrunoffPer'//trim(scratch)//trim(suffix)//grxt !4/26/2010 replaced "'.txt'" w/ "grxt"
	    outfil=trim(folder)//trim(rofil)
    	    call ssvgrd(ro,imax,pf1,nrow,ncol,u(3),test1,param,u(19),&
     &	    outfil,ti,header)
! save infiltration rate for grids that have non-zero runoff (changed unit number, from 6 to 8, 12/6/2010, RLB)
	    if (outp(6)) then
	      irfil='TRinfilratPer'//trim(scratch)//trim(suffix)//grxt
	      outfil=trim(folder)//trim(irfil)
   	      call ssvgrd(ir,imax,pf1,nrow,ncol,u(8),test1,param,u(19),&
     &	      outfil,ti,header)
	    end if
	   end if
! RLB--27 June 2011	  
	  else
	    write(scratch,'(i6)') j
	    scratch=adjustl(scratch)
	    if(rodoc) write(u(19),*) 'Runoff was zero at all grid cells for period ', scratch 
	    if(.not. rodoc) write(u(19),*) 'Runoff grid was not saved for period ', scratch 
	  end if  
   	end do
	else
! skip runoff routing if any of the input files are missing	
!  compute normalized infiltration intensities, Itransient/Ks, rik,
  	  write (*,*) 'Skipped runoff-routing computations;'
  	  write (*,*) 'Runoff routing input data did not exist.'
  	  write(u(19),*) ''
  	  write (u(19),*) 'Skipped runoff-routing computations;'
  	  write (u(19),*) 'Runoff routing input data did not exist.'
  	  write(u(19),*) '*******************************************'
! Added 2/2/2011 RLB
  	  nwf=1 
	  allocate (dsc(nwf),wf(nwf))
	do j=1,nper
!  read precipitation intensity, I
	 if (cri(j).lt.0) then
	  call srdgrd(grd,col,ncol,nrow,celsiz,nodat,ri,&
     &	  pf1,sctr,imax,temp,u(13),rifil(j),param,header,u(19))
    	  write(u(19),*) 'Precipitation intensity grid ',j
     	  write(u(19),*) trim(rifil(j)),sctr,' data cells'
     	if(sctr/=imx1) write (u(19),*) 'Grid mismatch ',trim(rifil(j))
  	 else
	  do  i=1,imx1
	   ri(i)=cri(j)
  	  end do
  	 end if
! compute infiltration at each cell for each time period  
	  do i=1,imx1
  	    if (ks(zo(i)).lt.ri(i)) then
	      ir(i)=ks(zo(i))
	      rik(i+(j-1)*imax)=1.d0
	    else
	      ir(i)=ri(i)
  	      rik(i+(j-1)*imax)=ri(i)/ks(zo(i))
	    end if
  	  end do
! save actual infiltration rates to files	  
! write only grids that are at or before the current time.
	  if (outp(6)) then
	   if(capt(j)<=t) then
	    ti=tiny(param(1))
	    write(scratch,'(i6)') j
	    scratch=adjustl(scratch)
	      irfil='TRinfilratPer'//trim(scratch)//trim(suffix)//grxt
	      outfil=trim(folder)//trim(irfil)
   	      call ssvgrd(ir,imax,pf1,nrow,ncol,u(8),test1,param,u(19),&
     &	      outfil,ti,header)
	    end if
	  end if
	end do
	end if ! end runouff routing and infiltration rates 
	return
  400	continue
  	write (*,*) '*** Error opening file in subroutine rnoff ***'
  	write (*,*) '--> ',infil
	write (*,*) 'Check file name and location'
  	write (u(19),*) '*** Error opening file in subroutine rnoff ***'
  	write (u(19),*) '--> ',infil
	write (u(19),*) 'Check file name and location'
  	write(*,*) 'Press RETURN to exit'
  	read*
	stop '400'
  430	continue
  	write (*,*) 'Attempt to read past end of file in subroutine rnoff'
	write (*,*) '--> ',infil
	write (*,*) 'Check file contents and value of Imax'
  	write (u(19),*) 'Attempt to read past end of file in subroutine rnoff'
	write (u(19),*) '--> ',infil
	write (u(19),*) 'Check file contents and value of Imax'
  	write(*,*) 'Press RETURN to exit'
  	read*
	stop '430'
	end subroutine rnoff
