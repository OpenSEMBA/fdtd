// Included inside FDTD_Solver in semba_fdtd.cpp.
// Mirrors the active time-loop structure from src/main/timestepping.F90.

void launch_simulation() {
    still_planewave_time = true;
    planewave_switched_off = false;
    const int report_step_interval = 500;
    const bool runtime_reporting_enabled = (numSteps >= report_step_interval);
    const double totalMcells = static_cast<double>(NX) * static_cast<double>(NY) *
                               static_cast<double>(NZ) / 1.0e6;
    const auto run_begin = std::chrono::steady_clock::now();
    auto last_report_time = run_begin;
    int last_report_step = 0;
#ifdef CompileWithMTLN
    if (isMpiRoot()) openMtlnObservation();
#endif
    for (n = 0; n <= numSteps; n++) {
        currentTime = n * dt;
        this->step();
        if (!runtime_reporting_enabled) {
            continue;
        }
        const bool is_report_step =
            (n > 0 && n % report_step_interval == 0) || (n == numSteps);
        if (!is_report_step) {
            continue;
        }
        const auto now = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration<double>(now - run_begin).count();
        double elapsed_chunk_s =
            std::chrono::duration<double>(now - last_report_time).count();
        if (elapsed_chunk_s <= 0.0) elapsed_chunk_s = 1.0e-12;

        const int reported_steps = std::max(1, n - last_report_step);
        const int total_steps = std::max(1, n);
        const double mcells_inst =
            (static_cast<double>(reported_steps) * totalMcells) / elapsed_chunk_s;
        const double mcells_avg =
            (static_cast<double>(total_steps) * totalMcells) /
            std::max(elapsed_s, 1.0e-12);
        const int next_info_step = std::min(numSteps, n + report_step_interval);

        std::cout << "Mins. since start  : "
                  << static_cast<int>(std::ceil(elapsed_s / 60.0)) << std::endl;
        std::cout << "Next info at step: " << next_info_step << std::endl;
        std::cout << std::uppercase << std::scientific << std::setprecision(8)
                  << "Total Mcells:    " << totalMcells << std::endl;
        std::cout << "Mcells/sec  : " << mcells_inst << " ("
                  << last_report_step << " to " << n << ")" << std::endl;
        std::cout << "Mcells/sec  : " << mcells_avg << " (" << 0 << " to " << n
                  << ")" << std::endl;
        std::cout << "__________________________________________" << std::endl;
        std::cout << std::defaultfloat;
        std::cout.flush();
        last_report_step = n;
        last_report_time = now;
    }
}

void launch() {
    launch_simulation();
}

void step() {
    const bool profile_enabled = runtimeProfile.enabled;
    const auto step_begin = profile_enabled
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};
    flushPlanewaveOff();

    // TODO: translate advanceAnisotropicE, advanceEDispersiveE, and multiport E.
    profileBlock(runtimeProfile.advanceE, [&]() { advanceE(); });
#ifdef CompileWithMTLN
    profileBlock(runtimeProfile.wiresE, [&]() { advanceMtlnE(); });
#endif
    profileBlock(runtimeProfile.wiresE, [&]() { advanceHollandWiresE(); });
    profileBlock(runtimeProfile.pmlE, [&]() { advancePmlE(); });
    profileBlock(runtimeProfile.sgbcE, [&]() { advanceSgbcE(); });
    profileBlock(runtimeProfile.lumpedE, [&]() { advanceLumpedE(); });
    applyPecE();
    profileBlock(runtimeProfile.planewaveE, [&]() { advancePlaneWaveE(); });
    profileBlock(runtimeProfile.nodalE, [&]() { advanceNodalE(); });

    profileBlock(runtimeProfile.mpiE, [&]() { flushMpiElectricFieldsOneAxis(); });

    // TODO: translate advanceAnisotropicH, advanceMDispersiveH, advanceNodalH, advanceWiresH, and multiport H.
    profileBlock(runtimeProfile.advanceH, [&]() { advanceH(); });
    profileBlock(runtimeProfile.pmlH, [&]() { advancePmlBodyH(); });
    profileBlock(runtimeProfile.pmlH, [&]() { advanceMagneticCpml(); });
    minusCloneMagneticPmc();
    cloneMagneticPeriodic();
    profileBlock(runtimeProfile.sgbcH, [&]() { advanceSgbcH(); });
    profileBlock(runtimeProfile.planewaveH, [&]() { advancePlaneWaveH(); });
    minusCloneMagneticPmc();
    cloneMagneticPeriodic();
    applyPecH();
    applyMurH();

    profileBlock(runtimeProfile.mpiH, [&]() { flushMpiMagneticFieldsOneAxis(); });

    profileBlock(runtimeProfile.sampling, [&]() {
        sampleProbes();
        sampleMovieProbes();
    });
    if (profile_enabled) {
        runtimeProfile.steps += 1;
        runtimeProfile.stepTotal += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - step_begin).count();
    }
}

void timestepping() {
    step();
}
