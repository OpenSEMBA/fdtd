integer function test_rotate_spacesteps() bind(C) result(err)
    use smbjson
    use rotate_testingTools
    use nfde_rotate_m
    
    integer(kind=4) :: flagRotation1=0, flagRotation2=0
    type(Desplazamiento), target :: des, oldDes
    type(MatrizMedios), target :: matriz, oldMatriz
    type(Parseador) :: parser

    call createDemoDesplazamiento(des)
    call createDemoMatrizMedios(matriz)
    oldDes = des
    oldMatriz = matriz
    parser%despl=>des
    parser%matriz=>matriz

    call rotate_generateSpaceSteps (parser, 1)
    err = 0
    !Check matriz total
    call expect_eq_int(err, oldMatriz%totalX, parser%matriz%totalZ)
    call expect_eq_int(err, oldMatriz%totalY, parser%matriz%totalX)
    call expect_eq_int(err, oldMatriz%totalZ, parser%matriz%totalY)

    !Check despl n
    call expect_eq_int(err, oldDes%nX, parser%despl%nZ)
    call expect_eq_int(err, oldDes%nY, parser%despl%nX)
    call expect_eq_int(err, oldDes%nZ, parser%despl%nY)

    !Check despl m1
    call expect_eq_int(err, oldDes%mx1, parser%despl%mz1)
    call expect_eq_int(err, oldDes%my1, parser%despl%mx1)
    call expect_eq_int(err, oldDes%mz1, parser%despl%my1)

    !Check despl m2
    call expect_eq_int(err, oldDes%mx2, parser%despl%mz2)
    call expect_eq_int(err, oldDes%my2, parser%despl%mx2)
    call expect_eq_int(err, oldDes%mz2, parser%despl%my2)

    !Check despl origin
    call expect_eq_real(err, oldDes%originx, parser%despl%originz)
    call expect_eq_real(err, oldDes%originy, parser%despl%originx)
    call expect_eq_real(err, oldDes%originz, parser%despl%originy)

    !Check despl des
    call expect_eq_real_vect(err, oldDes%desX, parser%despl%desZ)
    call expect_eq_real_vect(err, oldDes%desY, parser%despl%desX)
    call expect_eq_real_vect(err, oldDes%desZ, parser%despl%desY)

    if(err/=0) then
        write(*,*) "Failed rotation 1"
        write(*,*) err
        flagRotation1 = 1
        err = 0
    end if
    

    call createDemoDesplazamiento(des)
    call createDemoMatrizMedios(matriz)
    parser%despl=>des
    call rotate_generateSpaceSteps (parser, 2)

    !Check matriz total
    call expect_eq_int(err, oldMatriz%totalX, parser%matriz%totalY)
    call expect_eq_int(err, oldMatriz%totalY, parser%matriz%totalZ)
    call expect_eq_int(err, oldMatriz%totalZ, parser%matriz%totalX)

    !Check despl n
    call expect_eq_int(err, oldDes%nX, parser%despl%nY)
    call expect_eq_int(err, oldDes%nY, parser%despl%nZ)
    call expect_eq_int(err, oldDes%nZ, parser%despl%nX)

    !Check despl m1
    call expect_eq_int(err, oldDes%mx1, parser%despl%mY1)
    call expect_eq_int(err, oldDes%my1, parser%despl%mZ1)
    call expect_eq_int(err, oldDes%mz1, parser%despl%mX1)

    !Check despl m2
    call expect_eq_int(err, oldDes%mx2, parser%despl%mY2)
    call expect_eq_int(err, oldDes%my2, parser%despl%mZ2)
    call expect_eq_int(err, oldDes%mz2, parser%despl%mX2)

    !Check despl origin
    call expect_eq_real(err, oldDes%originx, parser%despl%originY)
    call expect_eq_real(err, oldDes%originy, parser%despl%originZ)
    call expect_eq_real(err, oldDes%originz, parser%despl%originX)

    !Check despl des
    call expect_eq_real_vect(err, oldDes%desX, parser%despl%desY)
    call expect_eq_real_vect(err, oldDes%desY, parser%despl%desZ)
    call expect_eq_real_vect(err, oldDes%desZ, parser%despl%desX)

    if(err/=0) then
        write(*,*) "Failed rotation 2"
        write(*,*) err
        flagRotation2 = 1
        err = 0
    end if

    if (flagRotation1 == 1 .or. flagRotation2 == 1) then 
        err = 1
    end if

end function

integer function test_rotate_generate_current_field_sources() bind(C) result(err)
end function

integer function test_rotate_generate_plane_waves() bind(C) result(err)
end function

integer function test_rotate_generate_box_sources() bind(C) result(err)
end function

integer function test_rotate_generate_fronteras() bind(C) result(err)
end function

integer function test_rotate_generate_pecs() bind(C) result(err)
end function

integer function test_rotate_generate_pmcs() bind(C) result(err)
end function

integer function test_rotate_generate_nonMetals() bind(C) result(err)
end function

integer function test_rotate_generate_anisotropics() bind(C) result(err)
end function

integer function test_rotate_generate_thin_wires() bind(C) result(err)
end function

integer function test_rotate_generate_slanted_wires() bind(C) result(err)
end function

integer function test_rotate_generate_thin_slots() bind(C) result(err)
end function

integer function test_rotate_generate_lossy_thin_surface() bind(C) result(err)
end function

integer function test_rotate_generate_fdms() bind(C) result(err)
end function

integer function test_rotate_generate_sondas() bind(C) result(err)
end function

integer function test_rotate_generate_mas_sondas() bind(C) result(err)
end function

integer function test_rotate_generate_bloque_probes() bind(C) result(err)
end function


integer function test_rotate_generate_volumic_probes() bind(C) result(err)
end function

integer function test_rotate_mpi() bind(c) result(err)
    use NFDETypes
    use smbjson
    use rotate_testingTools
    use nfde_rotate_m
    INTEGER (KIND=4) ::  mpidir    
    TYPE (coords) :: coordinates, oldCoordinates
    err = 0

    coordinates%Xi = 4
    coordinates%Xe = 5
    coordinates%Yi = 6
    coordinates%Ye = 7
    coordinates%Zi = 8
    coordinates%Ze = 9
    coordinates%Xtrancos = 1
    coordinates%Ytrancos = 2
    coordinates%Ztrancos = 3
    coordinates%Or = 2 

    oldCoordinates = coordinates
    mpidir = 1
    call rotatempi(mpidir, coordinates)



    call expect_eq_int(err, oldCoordinates%Xi, coordinates%Zi)
    call expect_eq_int(err, oldCoordinates%Xe, coordinates%Ze)
    call expect_eq_int(err, oldCoordinates%Yi, coordinates%Xi)
    call expect_eq_int(err, oldCoordinates%Ye, coordinates%Xe)
    call expect_eq_int(err, oldCoordinates%Zi, coordinates%Yi)
    call expect_eq_int(err, oldCoordinates%Ze, coordinates%Ye)
    call expect_eq_int(err, oldCoordinates%Xtrancos, coordinates%Ztrancos)
    call expect_eq_int(err, oldCoordinates%Ytrancos, coordinates%Xtrancos)
    call expect_eq_int(err, oldCoordinates%Ztrancos, coordinates%Ytrancos)
    call expect_eq_int(err, coordinates%Or, 1)

    if(err/=0) then
        write(*,*) "Failed rotation 1"
        write(*,*) err
        flagRotation1 = 1
        err = 0
    end if

    coordinates = oldCoordinates
    mpidir = 2
    call rotatempi(mpidir, coordinates)

    call expect_eq_int(err, oldCoordinates%Xi, coordinates%Yi)
    call expect_eq_int(err, oldCoordinates%Xe, coordinates%Ye)
    call expect_eq_int(err, oldCoordinates%Yi, coordinates%Zi)
    call expect_eq_int(err, oldCoordinates%Ye, coordinates%Ze)
    call expect_eq_int(err, oldCoordinates%Zi, coordinates%Xi)
    call expect_eq_int(err, oldCoordinates%Ze, coordinates%Xe)
    call expect_eq_int(err, oldCoordinates%Xtrancos, coordinates%Ytrancos)
    call expect_eq_int(err, oldCoordinates%Ytrancos, coordinates%Ztrancos)
    call expect_eq_int(err, oldCoordinates%Ztrancos, coordinates%Xtrancos)
    call expect_eq_int(err, coordinates%Or, 3)

    if(err/=0) then
        write(*,*) "Failed rotation 2"
        write(*,*) err
        flagRotation1 = 1
        err = 0
    end if

    if (flagRotation1 == 1 .or. flagRotation2 == 1) then 
        err = 1
    end if

end function


integer function test_rotate_mpi_scaled() bind(c) result(err)
    use NFDETypes
    use smbjson
    use rotate_testingTools
    use nfde_rotate_m
    INTEGER (KIND=4) ::  mpidir    
    TYPE (coords_SCALED) :: coordinates, oldCoordinates
    err = 0

    coordinates%Xi = 4
    coordinates%Xe = 5
    coordinates%Xc = 1.0
    coordinates%Yi = 6
    coordinates%Ye = 7
    coordinates%Yc = 2.0
    coordinates%Zi = 8
    coordinates%Ze = 9
    coordinates%Zc = 3.0
    coordinates%Or = 2 

    oldCoordinates = coordinates
    mpidir = 1
    call rotatempi_scaled(mpidir, coordinates)

    call expect_eq_int(err, oldCoordinates%Xi, coordinates%Zi)
    call expect_eq_int(err, oldCoordinates%Xe, coordinates%Ze)
    call expect_eq_real(err, oldCoordinates%Xc, coordinates%Zc)
    call expect_eq_int(err, oldCoordinates%Yi, coordinates%Xi)
    call expect_eq_int(err, oldCoordinates%Ye, coordinates%Xe)
    call expect_eq_real(err, oldCoordinates%Yc, coordinates%Xc)
    call expect_eq_int(err, oldCoordinates%Zi, coordinates%Yi)
    call expect_eq_int(err, oldCoordinates%Ze, coordinates%Ye)
    call expect_eq_real(err, oldCoordinates%Zc, coordinates%Yc)
    call expect_eq_int(err, coordinates%Or, 1)

    if(err/=0) then
        write(*,*) "Failed rotation 1"
        write(*,*) err
        flagRotation1 = 1
        err = 0
    end if

    coordinates = oldCoordinates
    mpidir = 2
    call rotatempi_scaled(mpidir, coordinates)

    call expect_eq_int(err, oldCoordinates%Xi, coordinates%Yi)
    call expect_eq_int(err, oldCoordinates%Xe, coordinates%Ye)
    call expect_eq_real(err, oldCoordinates%Xc, coordinates%Yc)
    call expect_eq_int(err, oldCoordinates%Yi, coordinates%Zi)
    call expect_eq_int(err, oldCoordinates%Ye, coordinates%Ze)
    call expect_eq_real(err, oldCoordinates%Yc, coordinates%Zc)
    call expect_eq_int(err, oldCoordinates%Zi, coordinates%Xi)
    call expect_eq_int(err, oldCoordinates%Ze, coordinates%Xe)
    call expect_eq_real(err, oldCoordinates%Zc, coordinates%Xc)
    call expect_eq_int(err, coordinates%Or, 3)

    if(err/=0) then
        write(*,*) "Failed rotation 2"
        write(*,*) err
        flagRotation1 = 1
        err = 0
    end if

    if (flagRotation1 == 1 .or. flagRotation2 == 1) then 
        err = 1
    end if

end function