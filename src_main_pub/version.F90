module  version 
character (128), parameter :: compdate= __DATE__ // " " // __TIME__  
character(128), parameter :: dataversion=&
'SEMBA-FDTD Version v1.0. Compiled on: ' // trim(adjustl(compdate))  
end module version
