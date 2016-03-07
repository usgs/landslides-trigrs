        subroutine dzero_brac(nstp,x,zroctr,zrptr)
! Procedure to bracket locations of zero in a list of values using change of sign between values
! Rex L Baum, USGS, 2 March 2012
        implicit none
	integer, parameter:: double=kind(1d0)
        integer:: n,nstp,zroctr,zrptr(nstp+1)
        real (double):: x(nstp+1)
        zroctr=0; zrptr=0
        do n=1,nstp+1 ! Check first for zero values in list
          if(x(n)==0.d0) then
            zroctr=zroctr+1
            zrptr(n)=2
          end if
        end do
        do n=1,nstp  ! Now check for zero values between points in list
          if(x(n)*x(n+1)<0.d0) then
            zroctr=zroctr+1
            zrptr(n)=1
          else
            cycle
          end if
        end do
        end subroutine dzero_brac