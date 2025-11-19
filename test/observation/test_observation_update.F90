integer function test_update_observation() bind(C) result(err)
    use observation_testingTools
    use OBSERVA
    use FDETYPES

    type(nf2ff_t) :: facesNF2FF
    type(sim_control_t) :: control

    type(taglist_t) :: tagNumbers
    type(media_matrices_t) :: media

    type(bounds_t)  ::  bounds
    type(limit_t), dimension(1:6) :: SINPML_fullsize
    
    logical :: ThereAreObservation, ThereAreWires, ThereAreFarFields
    integer :: initialtimestep, lastexecutedtime
    REAL(KIND=RKIND) :: eps00, mu00

    call initialize_nf2ff(facesNF2FF, &
        .FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .TRUE.)

    call initialize_control_flags(control,
        0, size, 3, 1077, &
        'frequencySlices', 'holland', &
        .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
        facesNF2FF)

    call initialize_tag_numbers(tagNumbers)

    call initialize_media(media)

    call initialize_bounds(bounds)

    call initialize_sinpml_fullsize(SINPML_fullsize)

    ThereAreObservation = .FALSE.
    ThereAreWires = .FALSE.
    ThereAreFarFields = .FALSE.

    eps00
    mu00

    call InitObservation(this%sgg,media,tagNumbers, &
                                ThereAreObservation, ThereAreWires, ThereAreFarFields, initialtimestep, lastexecutedtime, &
                                this%sinPML_fullsize,this%eps0, this%mu0, bounds, control)
    
    call UpdateObservation()

end function test_update_observation