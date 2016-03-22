!-----------------------------------------------------------
!
! partial_p by M. Alvioli, November 2014  
!
!-----------------------------------------------------------
subroutine partial_p(imax)
  use mpi
  use partial_arrays_p
  implicit none
  integer Nmin,Next,k,i,imax
  integer myrank,isize,ierr
  !
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  allocate(isct(0:isize-1),idsp(0:isize-1))
  !
  Nmin=imax/isize
  Next=mod(imax,isize)
  k=0
  DO i=0,isize-1
     if(i.lt.Next) then
        isct(i)=Nmin+1
     else
        isct(i)=Nmin
     endif
     idsp(i)=k
     k=k+isct(i)
  enddo
  irct=isct(myrank)
  !
!  print *
!  print *,myrank,isct(myrank),idsp(myrank)
  !
  return
end subroutine partial_p 
!-----------------------------------------------------------
