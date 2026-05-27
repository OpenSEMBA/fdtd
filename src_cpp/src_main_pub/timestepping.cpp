// Included inside FDTD_Solver in semba_fdtd.cpp.
// Mirrors the active time-loop structure from src_main_pub/timestepping.F90.

void launch_simulation() {
    still_planewave_time = true;
    planewave_switched_off = false;
    std::cout << "Running FDTD: " << numSteps << " steps..." << std::endl;
#ifdef CompileWithMTLN
    if (isMpiRoot()) openMtlnObservation();
#endif
    for (n = 0; n <= numSteps; n++) {
        currentTime = n * dt;
        this->step();
        if (n % 500 == 0 || n == numSteps)
            std::cout << "  Step " << n << "/" << numSteps << " (t=" << currentTime << "s)" << std::endl;
    }
    std::cout << "Simulation complete." << std::endl;
}

void launch() {
    launch_simulation();
}

void step() {
    flushPlanewaveOff();

    // TODO: translate advanceAnisotropicE, advanceEDispersiveE, and multiport E.
    advanceE();
#ifdef CompileWithMTLN
    advanceMtlnE();
#endif
    advanceHollandWiresE();
    advancePmlE();
    advanceSgbcE();
    advanceLumpedE();
    applyPecE();
    advancePlaneWaveE();
    applyPecE();
    advanceNodalE();

    flushMpiElectricFieldsOneAxis();

    // TODO: translate advanceAnisotropicH, advanceMDispersiveH, advanceNodalH, advanceWiresH, and multiport H.
    advanceH();
    advancePmlBodyH();
    advanceMagneticCpml();
    minusCloneMagneticPmc();
    cloneMagneticPeriodic();
    advanceSgbcH();
    advancePlaneWaveH();
    minusCloneMagneticPmc();
    cloneMagneticPeriodic();
    applyPecH();
    applyMurH();

    flushMpiMagneticFieldsOneAxis();

    sampleProbes();
    sampleMovieProbes();
}

void timestepping() {
    step();
}
