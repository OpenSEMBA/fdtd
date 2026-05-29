#include "observation_preprocess.h"
#include "nfde_types.h"

namespace Observa_m {

using NFDETypes_m::mapvtk;

void preprocess_observation_full(Obses_t_full& observation, output_t& privateOutput,
                                 const std::vector<RKIND_tiempo>& time,
                                 int finaltimestep, RKIND_tiempo dt, bool saveall) {
    observation.done = false;
    observation.flushed = false;
    observation.begun = false;

    observation.TimeStep = std::max(observation.TimeStep, dt);

    const RKIND span = observation.FinalTime - observation.InitialTime;
    const RKIND minStep = std::min(dt, observation.TimeStep);
    if (10.0 * span / minStep >= static_cast<RKIND>(std::numeric_limits<int32_t>::max())) {
        observation.FinalTime = observation.InitialTime +
            minStep * static_cast<RKIND>(std::numeric_limits<int32_t>::max()) / 10.0;
    }

    if (observation.InitialTime < observation.TimeStep) {
        observation.InitialTime = 0.0;
    }

    if (observation.TimeStep > (observation.FinalTime - observation.InitialTime)) {
        if (!observation.P.empty() && observation.P[0].what == mapvtk) {
            observation.FinalTime = 0.0;
            observation.InitialTime = 0.0;
        } else {
            observation.FinalTime = observation.InitialTime + observation.TimeStep;
        }
    }

    observation.FreqStep = std::min(observation.FreqStep, 2.0 / dt);
    if ((observation.FreqStep > observation.FinalFreq - observation.InitialFreq) ||
        observation.FreqStep == 0.0) {
        observation.FreqStep = observation.FinalFreq - observation.InitialFreq;
        observation.FinalFreq = observation.InitialFreq + observation.FreqStep;
    }

    if (!observation.Volumic) {
        observation.Saveall = observation.Saveall || saveall;
        privateOutput.SaveAll = observation.Saveall;
    } else {
        privateOutput.SaveAll = false;
        observation.Saveall = false;
    }

    if (!observation.P.empty() && observation.P[0].what == mapvtk) {
        privateOutput.SaveAll = false;
        observation.Saveall = false;
    }

    if (observation.Saveall) {
        privateOutput.Trancos = 1;
        observation.InitialTime = 0.0;
        if (finaltimestep + 2 < static_cast<int>(time.size())) {
            observation.FinalTime = time[static_cast<size_t>(finaltimestep + 2)];
        }
    } else {
        privateOutput.Trancos = std::max(1, static_cast<int>(observation.TimeStep / dt));
        observation.InitialTime = std::max(0.0, observation.InitialTime);
        if (finaltimestep + 2 < static_cast<int>(time.size())) {
            observation.FinalTime =
                std::min(time[static_cast<size_t>(finaltimestep + 2)], observation.FinalTime);
        }
        if (observation.FinalTime < observation.InitialTime) {
            observation.FinalTime = observation.InitialTime;
        }
    }
}

} // namespace Observa_m
