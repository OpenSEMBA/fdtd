integer function test_rotate_spacesteps() bind(C) result(err)
    use smbjson
    use rotate_testingTools
    use nfde_rotate_m

    type(Desplazamiento) :: des, oldDes
    type(Parseador) :: parser
    call createDemoDesplazamiento(des)
    oldDes = des
    parser%despl=des

    call rotate_generateSpaceSteps (parser, 1)

    err = 0
    call expect_eq_real_vect(err, oldDes%desY, parser%despl%desX)

end function
