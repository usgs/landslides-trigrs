	subroutine mpfldr(rc,dir,nodata)
c ! Maps ESRI arc/grid flow directions to trigrs flow directions
c ! by Rex L. Baum, USGS, 22 April 2002
	implicit none
	integer i,rc,dir(rc),adr(255),nodata
	do 10, i=1,255
c ! Cells having undefined flow directions point to themselves
	  adr(i)=5
   10	continue
c ! Define correspondence between flow directions
   	adr(32)=1
   	adr(64)=2
   	adr(128)=3
   	adr(16)=4
   	adr(1)=6
   	adr(8)=7
   	adr(4)=8
   	adr(2)=9
	do 20,i=1,rc
	if(dir(i).ne.nodata) dir(i)=adr(dir(i))
   20	continue
   	return
	end  
	