//#include "plot_histogram.hpp"
#include "plot_simulation.h"
#include "track_resolution.h"
#include "plot_MobilityModel.h"
#include "plot_ComparingVoltage.h"
#include "plot_DriftTime.h"
#include "plot_inPixel.h"
#include <time.h>

// plot_histogram::plot_histogram() {
//     std::cout << "plot_histogram object is created" << std::endl;
// }

// plot_simulation::plot_simulation() {
//     std::cout << "plot_simulation object is created" << std::endl;
// }

int main() {
    clock_t start = clock();
    
    gErrorIgnoreLevel = kWarning;

    plot_simulation plot;
    track_resolution track;
    plot_MobilityModel plot_Mobility;
    plot_ComparingVoltage plot_voltage;
    plot_DriftTime plot_driftTime;
    plot_inPixel plot_inPixel;
    // plot_Mobility.run();
    // plot_voltage.voltage_run();
    plot_driftTime.run_driftTime();
    plot_inPixel.run_inPixel();

    LOG_STATUS.source("ce65.cpp/main") << "Main process  is complete.";

    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
    std::cout << "Time: " << time << "ms" << std::endl;
}