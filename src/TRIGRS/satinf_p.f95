subroutine satinf_p(imx1,ulog,u1,nccs)
  ! 7/20/2006 Rex L. Baum, USGS, Latest revision 29 Jan 2013 
  !  computations for time series in the fully saturated case.
  !  calls ivestp() for infinite depth solution.	
  !  
  ! MPI calls & MPI-related changes by M. Alvioli, November 2014
  !
  use grids; use input_vars 
  use model_vars
  use input_file_defs
  use mpi
  use modules_p
  use partial_arrays_p
  implicit none
  integer::i,j,jf,k,ulog,u1,imx1,nccs
  !	integer::nmn1,nmin1,nmax3,nmax0,nmns,nmxs
  !	logical:: lcvs 
  real:: qbij(nts+1)
  real (double)::rf(nzs+1),finf
  integer myrank,isize,ierr
  integer nccs_0
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  !	nmax3=0;nmax0=0
  !	nmn1=nmax+1;nmin1=nmax+1
  if(myrank.eq.0) then
     write(ulog,*) 'Starting computations'
     write(ulog,*) 'for infinite-depth saturated zone'
     write(*,*) 'Starting computations'
     write(*,*) 'for infinite-depth saturated zone'
     write(*,*) 'Cells completed: '
  endif
  ! loop over all grid cells
  finf=10.
  grid_loop: do i=idsp(myrank)+1,idsp(myrank)+isct(myrank)
     if(i.gt.imx1) goto 999
     if (mod(i-1,2000)==0.and.myrank.eq.0) write (*,'(2x,f5.1,a1,2x,a1)',advance="no")&
          &float(i-idsp(myrank))/float(irct)*100,"%",char(13) ! cells completed
     if(slo(i)<slomin .or. slo(i)>slomax .or. zmax(i)<=0.0001) then ! default values for gently or steeply sloping cells 
        do jf=1,nout
           fsmin(i+(jf-1)*isct(myrank)-idsp(myrank))=finf+1.
           zfmin(i+(jf-1)*isct(myrank)-idsp(myrank))=zmax(i)
           pmin(i+(jf-1)*isct(myrank)-idsp(myrank))=0.
        end do
        cycle
     end if
     !	    lcvs=.true.
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
     call ivestp_p(u1,qbij,ulog,i,rf) ! Revised 28 Jan 2013, RLB 
     !                & ulog,i,rf,nccs,lcvs,nmxs)
     !                nmns=nmn
999  continue
  end do grid_loop
  if(myrank.eq.0) then
     if(myrank.eq.0) write(*,*)
     write (*,*) imx1, ' cells completed' 
     write (ulog,*) imx1, ' cells completed' 
  endif
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myrank.eq.0) write(*,*) " (MPI) collecting results .. "
  if(outp(3)) then
     do j=1,nout
        Call MPI_GATHERV(fsmin(1+irct*(j-1):irct*j),irct,MPI_FLOAT,&
             &fsmin_0(1+(j-1)*imx1),isct,idsp,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
     enddo
  endif
  if(outp(4)) then
     do j=1,nout
        Call MPI_GATHERV(zfmin(1+irct*(j-1):irct*j),irct,MPI_FLOAT,&
             &zfmin_0(1+(j-1)*imx1),isct,idsp,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
     enddo
  endif
  if(outp(5)) then
     do j=1,nout
        Call MPI_GATHERV(pmin(1+irct*(j-1):irct*j),irct,MPI_FLOAT,&
             &pmin_0(1+(j-1)*imx1),isct,idsp,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
     enddo
  endif
  return
end subroutine satinf_p
