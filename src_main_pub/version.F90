module  version 
    character(128), parameter :: compilation_date= __DATE__ // " " // __TIME__  
#ifdef CompileWithDebug
    character(128), parameter :: compilation_mode='Debug build'
#else
    character(128), parameter :: compilation_mode='Release build' 
#endif
    character(128), parameter :: dataversion = &
        'SEMBA-FDTD. '// trim(adjustl(compilation_mode)) // &
            ', compiled on: ' // trim(adjustl(compilation_date))  
end module version
