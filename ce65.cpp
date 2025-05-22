//#include "plot_histogram.hpp"
#include "plot_simulation.hpp"
#include "track_resolution.hpp"
#include "plot_MobilityModel.hpp"
#include "plot_ComparingVoltage.hpp"
#include <time.h>

// plot_histogram::plot_histogram() {
//     std::cout << "plot_histogram object is created" << std::endl;
// }

// plot_simulation::plot_simulation() {
//     std::cout << "plot_simulation object is created" << std::endl;
// }

int main() {
    clock_t start = clock();

    plot_simulation plot;
    track_resolution track;
    plot_MobilityModel plot_Mobility;
    plot_ComparingVoltage plot_voltage;
    plot_Mobility.run();
    plot_voltage.voltage_run();

    std::cout << "Main function is executed" << std::endl;

    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
    std::cout << "Time: " << time << "ms" << std::endl;
}