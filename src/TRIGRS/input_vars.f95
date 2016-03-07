! By Rex L. Baum, 1 April 2004
! 19Aug2009 RLB added bkgrof
module input_vars
	logical:: ans,outp(8),rodoc,lskip,lany,llus,lps0,unsat0,bkgrof
	logical:: lasc,lpge0 ! Added 4/14-15/2010
	logical:: igcapf ! added 2/15/2012 RLB
	logical,allocatable:: unsat(:), igcap(:)  ! added igcap(:) 2/15/2012 RLB
 	integer:: imax,row,col,nwf,tx,nmax
	integer:: flag,nper,spcg ! Added spcg 2/14/2012 RLB
  	integer:: nzs,mmax
	integer:: nzon,nout
	integer,allocatable:: ksav(:), uijz(:) !added uijz, 12/7/2010 RLB
	real:: uww,zmin,t,dep,czmax,crizero,slomin,slomax,deepz
	real,allocatable:: ths(:),thr(:),alp(:),dif(:),c(:),phi(:)
	real,allocatable:: ks(:),uws(:),capt(:),cri(:),tsav(:)
	character (len=5):: flowdir, el_or_dep
	character (len=4):: deepwat
end module input_vars
