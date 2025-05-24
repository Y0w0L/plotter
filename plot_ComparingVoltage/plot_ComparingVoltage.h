#ifndef PLOT_COMPARINGVOLTAGE_H
#define PLOT_COMPARINGVOLTAGE_H

#include "plot_histogram/plot_histogram.h"
#include "plot_MobilityModel/plot_MobilityModel.h"

class plot_ComparingVoltage : public plot_MobilityModel {
    public:
        plot_ComparingVoltage();
        void get_voltageClsize(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage);
        void get_voltageCharge(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage);
        void get_voltageResidual(std::vector<TFile*> input_ROOTFile, TFile* output, std::vector<std::string> voltage);
        void voltage_run();
};

#endif