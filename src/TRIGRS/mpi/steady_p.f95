!
! MPI calls & MPI-related changes by M. Alvioli, November 2014  
!
subroutine steady_p(sumex,ulog,imx1)
  use grids; use input_vars
  use mpi
  implicit none
  integer:: i,acnt,ulog,imx1
  integer,parameter:: double=kind(1d0)
  real (double):: sumex,b,rslo
  integer myrank,isize,ierr
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  !  By Rex L. Baum, 1 April 2004
  !  Compute initial estimate of Isteady/Ks, rikzero(), and test values.
  !  Isteady must be < Ks.  If cos(slo) *cos(slo)<(Isteady/Ks), then an 
  !  inverted water table results.
  if(myrank.eq.0) write(*,*) 'Testing and adjusting steady infiltration rates'
  acnt=0
  do i=1,imx1
     if(ks(zo(i))==0.) then ! prevent division by zero errors
        rikzero(i)=1.
     else
        rikzero(i)=rizero(i)/ks(zo(i))
     end if
  end do
  sumex=0
  do i=1,imx1
     rslo=slo(i)
     b=cos(rslo)*cos(rslo)
     if(rikzero(i)>=b) then
        rikzero(i)=cos(rslo) 
        acnt=acnt+1
        if(myrank.eq.0) then
           write(ulog,*) 'Adjusted steady infiltration rate, cell ',i
           write(ulog,*) 'Corrected rizero/Ks = ',rikzero(i)
           write(ulog,*) 'Original rizero, ks',rizero(i),ks(zo(i))
           write(ulog,*) 'Set pore presssures to zero'
        endif
     end if
     if (depth(i)==0 .and. rizero(i)<0) then
        sumex=sumex-rizero(i)
     end if
  end do
  if(myrank.eq.0) then
     write(*,*) 'Adjusted steady infiltration rate at ',acnt,' cells'
     write(ulog,*) 'Adjusted steady infiltration rate at ',acnt,' cells'
  endif
  return
end subroutine steady_p
