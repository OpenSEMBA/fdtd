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
    function create_tag_list(sggAlloc) result(r)
        type(XYZlimit_t), dimension(6), intent(in) :: sggAlloc
        type(taglist_t) :: r


        allocate (r%edge%x(sggAlloc(iEx)%XI:sggAlloc(iEx)%XE, sggAlloc(iEx)%YI:sggAlloc(iEx)%YE, sggAlloc(iEx)%ZI:sggAlloc(iEx)%ZE))
        allocate (r%edge%y(sggAlloc(iEy)%XI:sggAlloc(iEy)%XE, sggAlloc(iEy)%YI:sggAlloc(iEy)%YE, sggAlloc(iEy)%ZI:sggAlloc(iEy)%ZE))
        allocate (r%edge%z(sggAlloc(iEz)%XI:sggAlloc(iEz)%XE, sggAlloc(iEz)%YI:sggAlloc(iEz)%YE, sggAlloc(iEz)%ZI:sggAlloc(iEz)%ZE))
        allocate (r%face%x(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHx)%YI:sggAlloc(iHx)%YE, sggAlloc(iHx)%ZI:sggAlloc(iHx)%ZE))
        allocate (r%face%y(sggAlloc(iHy)%XI:sggAlloc(iHy)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHy)%ZI:sggAlloc(iHy)%ZE))
        allocate (r%face%z(sggAlloc(iHz)%XI:sggAlloc(iHz)%XE, sggAlloc(iHz)%YI:sggAlloc(iHz)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))


        r%edge%x(:,:,:) = 0
        r%edge%y(:,:,:) = 0
        r%edge%z(:,:,:) = 0
        r%face%x(:,:,:) = 0
        r%face%y(:,:,:) = 0
        r%face%z(:,:,:) = 0
    end function create_tag_list

    function create_media(sggAlloc) result(r)
        type(XYZlimit_t), dimension(6), intent(in) :: sggAlloc
        type(media_matrices_t) :: r

        allocate (r%sggMtag(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))
        allocate (r%sggMiNo(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))

        allocate (r%sggMiEx(sggAlloc(iEx)%XI:sggAlloc(iEx)%XE, sggAlloc(iEx)%YI:sggAlloc(iEx)%YE, sggAlloc(iEx)%ZI:sggAlloc(iEx)%ZE))
        allocate (r%sggMiEy(sggAlloc(iEy)%XI:sggAlloc(iEy)%XE, sggAlloc(iEy)%YI:sggAlloc(iEy)%YE, sggAlloc(iEy)%ZI:sggAlloc(iEy)%ZE))
        allocate (r%sggMiEz(sggAlloc(iEz)%XI:sggAlloc(iEz)%XE, sggAlloc(iEz)%YI:sggAlloc(iEz)%YE, sggAlloc(iEz)%ZI:sggAlloc(iEz)%ZE))

        allocate (r%sggMiHx(sggAlloc(iHx)%XI:sggAlloc(iHx)%XE, sggAlloc(iHx)%YI:sggAlloc(iHx)%YE, sggAlloc(iHx)%ZI:sggAlloc(iHx)%ZE))
        allocate (r%sggMiHy(sggAlloc(iHy)%XI:sggAlloc(iHy)%XE, sggAlloc(iHy)%YI:sggAlloc(iHy)%YE, sggAlloc(iHy)%ZI:sggAlloc(iHy)%ZE))
        allocate (r%sggMiHz(sggAlloc(iHz)%XI:sggAlloc(iHz)%XE, sggAlloc(iHz)%YI:sggAlloc(iHz)%YE, sggAlloc(iHz)%ZI:sggAlloc(iHz)%ZE))

        r%sggMtag (:, :, :) = 0 
        r%sggMiNo (:, :, :) = 1
        r%sggMiEx (:, :, :) = 1
        r%sggMiEy (:, :, :) = 1
        r%sggMiEz (:, :, :) = 1
        r%sggMiHx (:, :, :) = 1
        r%sggMiHy (:, :, :) = 1
        r%sggMiHz (:, :, :) = 1
    end function create_media

end module FDETYPES_TOOLS