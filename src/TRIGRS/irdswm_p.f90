!  
!  Reads a list of integer values that represent a 2-D matrix into 
!  two 1-D arrays.  One array is a pointer that tracks the starting
!  location of each row in the matrix.  This scheme is advantageous 
!  for storing sparse "sawtooth" arrays. Such arrays have rows of 
!  variable length and all of the non-zero values are at the left 
!  ends of the rows. 
!  by Rex L. Baum, spring 2001 
!  
! Implementation:
!  	open (u,file=infil,status='old')
!  	call irdswm(len,jmax,u,test,x,u2,ctr)
!  	close(u)
!
! MPI calls & MPI-related changes by M. Alvioli, November 2014
!  
subroutine irdswm_p(len,jmax,u,test,x,ctr,ulog)
  use mpi
  implicit none
  integer j,k,jmax,u,val,test,len,x(len),ulog
  integer ctr(jmax+1)
  integer myrank,isize,ierr,jump
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
  j=0
  ctr(1)=0
  k=0
  jump=0
  if(myrank.eq.0) then
100  continue
     read (u,*,end=110) val
     !  Check for marker value at end of a row
     if (val.eq.test) then
        read (u,*,end=110) j
        ctr(j)=k+1
        if (j.gt.jmax) then
           jump=1
           go to 150
        endif
        go to 100
     end if
     k=k+1
     x(k)=val
     go to 100
110  continue
     ctr(jmax+1)=k+1
  endif
150 continue
  Call MPI_BCAST(jump,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(jump.eq.1) goto 200
  Call MPI_BCAST(ctr,jmax+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call MPI_BCAST(x,len,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  return
200 continue
  if(myrank.eq.0) then
     write(*,*) 'Array bound <jmax> exceeded in subroutine irdswm'
     write(ulog,*) 'Array bound <jmax> exceeded in subroutine irdswm'
  endif
  call MPI_FINALIZE(ierr)
  stop '-200 in irdswm()'
end subroutine irdswm_p
