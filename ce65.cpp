//#include "plot_histogram.hpp"
#include "plot_simulation.h"
#include "track_resolution.h"
#include "plot_MobilityModel.h"
#include "plot_ComparingVoltage.h"
#include "plot_DriftTime.h"
#include "plot_inPixel.h"
#include "plot_LabTest.h"
#include "plot_ExperimentData.h"
#include "plot_BeamTest.h"
#include "plot_AllpixData.h"
#include <time.h>

// plot_histogram::plot_histogram() {
//     std::cout << "plot_histogram object is created" << std::endl;
// }

// plot_simulation::plot_simulation() {
//     std::cout << "plot_simulation object is created" << std::endl;
// }

int main(int argc, char* argv[]) {
    clock_t start = clock();
    
    gErrorIgnoreLevel = kWarning;

    plot_simulation plot;
    track_resolution track_resolution;
    plot_MobilityModel plot_Mobility;
    plot_ComparingVoltage plot_voltage;
    plot_DriftTime plot_driftTime;
    plot_inPixel plot_inPixel;
    plot_ExperimentData plot_ExperimentData;
    plot_LabTest plot_LabTest;
    plot_BeamTest plot_BeamTest;
    //allpix::plot_AllpixData plot_AllpixData;
    // plot_Mobility.run();
    // plot_voltage.voltage_run();
    // plot_driftTime.run_driftTime();
    //plot_inPixel.run_inPixel();
    //track_resolution.run();
    //plot_ExperimentData.run_inPixel();
    //plot_ExperimentData.run_Analysis();
    //plot_ExperimentData.run_NoiseScan();
    //plot_LabTest.run_LabTest();

    //
    // plot_BeamTest kek202412_plotter(DataSource::KEK202412);
    // kek202412_plotter.run_kek_plots();
    // plot_BeamTest sps202404_plotter(DataSource::SPS202404);
    // sps202404_plotter.run_sps_plots();

    try {
        plot_BeamTest.BeamTest_main(argc, argv);
    } catch (const std::exception& e) {
        LOG_ERROR.source("ce65.cpp/main") << "ERROR happens dourning processing. " << e.what() ;
        return 1;
    }
    //plot_AllpixData.run_analysis();

    LOG_STATUS.source("ce65.cpp/main") << "Main process finished";

    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
    std::cout << "Time: " << time << "ms" << std::endl;

    return 0;
}