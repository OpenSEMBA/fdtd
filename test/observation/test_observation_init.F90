integer function test_init_time_movie_observation() bind(C) result(err)
    use FDETYPES
    use FDETYPES_TOOLS
    use Observa
    use observation_testingTools

    type(SGGFDTDINFO) ::  sgg
    type(media_matrices_t) :: media
    type(taglist_t) :: tag_numbers
    logical :: ThereAreObservation, ThereAreWires, ThereAreFarFields
    integer :: initialtimestep
    real(kind=RKIND_tiempo) :: lastexecutedtime
    type(limit_t), dimension(1:6) :: SINPML_fullsize
    type(bounds_t)  ::  bounds
    type(sim_control_t) :: control
    real(kind=RKIND) :: eps = EPSILON_VACUUM, mu = MU_VACUUM

    type(output_t), pointer, dimension(:) :: output

    sgg = create_base_sgg()
    SINPML_fullsize = create_limit_t(0,4,0,4,0,4,3,3,3)

    call InitObservation(sgg, media, tag_numbers, &
        ThereAreObservation, ThereAreWires, ThereAreFarFields, initialtimestep, lastexecutedtime, &
        SINPML_fullsize, eps, mu, bounds, control)
    output => GetOutput()

    !assertion
end function test_init_time_movie_observation