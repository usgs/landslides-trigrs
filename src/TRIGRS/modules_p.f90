!-----------------------------------------------------------
!
! modules_p by M. Alvioli, November 2014  
!
!-----------------------------------------------------------
Module modules_p
  integer,allocatable:: nv_0(:),nvu_0(:)
  real,allocatable :: zfmin_0(:),fsmin_0(:),pmin_0(:)
  real,allocatable :: wtab_0(:)
  double precision,allocatable:: p3d_0(:,:),pzero3d_0(:,:),ptran3d_0(:,:),fs3d_0(:,:),th3d_0(:,:),dh3d_0(:),newdep3d_0(:)
End Module modules_p
!-----------------------------------------------------------
!
! addition March 2015
!
!-----------------------------------------------------------
Module partial_arrays_p
  integer :: irct
  allocatable :: isct(:),idsp(:)
End Module partial_arrays_p
!-----------------------------------------------------------
