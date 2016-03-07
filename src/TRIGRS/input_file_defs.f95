module input_file_defs
  	character (len=255):: heading,slofil,zonfil,zfil
  	character (len=255):: depfil
  	character (len=255):: rizerofil,nxtfil,ndxfil
	character (len=255):: wffil,dscfil,init,title
  	character (len=255):: elevfil ! added 4/21/2010
  	character (len=255), allocatable:: rifil(:)
 	character (len=224):: folder,elfoldr
	character (len=8):: suffix
end module input_file_defs

