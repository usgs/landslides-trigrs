! By Rex L. Baum, 1 April 2004
module grids
	integer,allocatable:: pf2(:),indx(:),nxt(:),nv(:),nvu(:)
	integer,allocatable:: dsctr(:),dsc(:),zo(:), itemp(:) ! Added itemp() 2/13/2013
	real,allocatable:: rikzero(:)
	real,allocatable::  rik(:),rik1(:),ri(:),rizero(:),pf1(:)
	real,allocatable:: temp(:),ro(:),wf(:),ir(:),tfg(:)
	real,allocatable:: zmax(:),slo(:),depth(:)
	real,allocatable:: zfmin(:),fsmin(:),pmin(:) 
	real,allocatable:: elev(:),wtab(:) ! added 4/21/2010, wtab() added 4/19/11
	character (len=4):: grxt ! Added grxt 4/26/2010	
end module grids