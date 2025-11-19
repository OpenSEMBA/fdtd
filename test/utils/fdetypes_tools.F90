module FDETYPES_TOOLS
    use FDETYPES
    contains
    function create_limit_t(XI,XE,YI,YE,ZI,ZE,NX,NY,NZ) result(r)
        type(limit_t) :: r
        integer (kind=4), intent(in) :: XI,XE,YI,YE,ZI,ZE,NX,NY,NZ
        r%XI = XI
        r%XE = XE
        r%YI = YI
        r%YE = YE
        r%ZI = ZI
        r%ZE = ZE
        r%NX = NX
        r%NY = NY
        r%NZ = NZ
    end function create_limit_t
end module FDETYPES_TOOLS