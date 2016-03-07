        subroutine svlist(u1)
	use input_file_defs; use input_vars
        use grids; use model_vars
        implicit none 
        integer::u1,i,j,n
        real (double):: ddg2rad,rslo
	real (double):: zns,zinc,z,b1
	pi=3.141592653589793
	ddg2rad=pi/180.D0
	if(flag==-1) then ! write simple list file with pressure head and factor of safety,
          grid_loop_1: do i=1,imax
	    rslo=slo(i)
	    time_loop_1: do n=1,nout
	      write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
	      &'cell ',i,rslo/ddg2rad,n,tsav(n) 
	      zns=float(nzs)
	      z=zmin
              zinc=(zmax(i)-zmin)/zns
	      Z_loop_1: do j=1,nzs+1
	        write(u1,'(6(g12.5,1x):)') z,p3d(i+(n-1)*imax,j),fs3d(i+(n-1)*imax,j)
	        z=z+zinc	  
	      end do Z_loop_1
            end do time_loop_1
          end do grid_loop_1
        else if (flag==-2) then ! write detailed list file, including initial, transient, and beta-line pressure head values.
          grid_loop_2: do i=1,imax
	    rslo=slo(i)
	    time_loop_2: do n=1,nout
	      write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
	      &'cell ',i,rslo/ddg2rad,n,tsav(n) 
	      zns=float(nzs)
	      z=zmin
              zinc=(zmax(i)-zmin)/zns
	      Z_loop_2: do j=1,nzs+1
	        b1=cos(rslo)
	        select case (flowdir) ! set value of beta (Iverson's beta line)
	          case ('slope')
	            beta=b1*b1
	          case ('hydro')
	            beta=1.d0
	          case default
	            beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
	        end select
	        if(abs(b1-rikzero(i))<1.e-6) beta=0.d0 
	        bline(j)=z*beta
                write(u1,'(6(g12.5,1x):)') z,p3d(i+(n-1)*imax,j),pzero3d(i,j),&
                &ptran3d(i+(n-1)*imax,j),bline(j),fs3d(i+(n-1)*imax,j)
	        z=z+zinc	  
	      end do Z_loop_2
            end do time_loop_2
            end do grid_loop_2
        else if(flag==-3) then ! Write list file with pressure head, factor of safety, and water content
          if (unsat0) then
            grid_loop_3a: do i=1,imax
	      rslo=slo(i)
	      time_loop_3a: do n=1,nout
	        write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
	        &'cell ',i,rslo/ddg2rad,n,tsav(n) 
	        zns=float(nzs)
	        z=zmin
                zinc=(zmax(i)-zmin)/zns
	        Z_loop_3a: do j=1,nzs+1
                  write(u1,'(6(g12.5,1x):)') z,p3d(i+(n-1)*imax,j),fs3d(i+(n-1)*imax,j),th3d(i+(n-1)*imax,j)
	          z=z+zinc	  
	        end do Z_loop_3a
              end do time_loop_3a
            end do grid_loop_3a
          else
            grid_loop_3b: do i=1,imax
	      rslo=slo(i)
	      time_loop_3b: do n=1,nout
	        write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
	        &'cell ',i,rslo/ddg2rad,n,tsav(n) 
	        zns=float(nzs)
	        z=zmin
                zinc=(zmax(i)-zmin)/zns
	        Z_loop_3b: do j=1,nzs+1
                  write(u1,'(6(g12.5,1x):)') z,p3d(i+(n-1)*imax,j),fs3d(i+(n-1)*imax,j)
	          z=z+zinc	  
	        end do Z_loop_3b
              end do time_loop_3b
            end do grid_loop_3b
          end if
        end if
        end subroutine svlist
