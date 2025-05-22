#ifndef PLOT_COMPARINGVOLTAGE_HPP
#define PLOT_COMPARINGVOLTAGE_HPP

#include "plot_histogram/plot_histogram.hpp"
#include "plot_MobilityModel/plot_MobilityModel.hpp"

class plot_ComparingVoltage : public plot_MobilityModel {
    public:
        plot_ComparingVoltage();
        void get_voltageClsize(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage);
        void get_voltageCharge(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage);
        void get_voltageResidual(std::vector<TFile*> input_ROOTFile, TFile* output, std::vector<std::string> voltage);
        void voltage_run();
};

#endif