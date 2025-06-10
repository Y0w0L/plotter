#ifndef PLOT_DRIFTTIME_H
#define PLOT_DRIFTTIME_H

#include "plot_histogram/plot_histogram.h"

class plot_DriftTime : public plot_histogram {
public:
    plot_DriftTime();
    static TCanvas* get_driftTime(std::vector<TFile*> inputROOTFile, TFile* outputROOTFile, const std::vector<std::string>& position_name);
    void run_driftTime();
};

#endif // PLOT_DRIFTTIME_H