#ifndef PLOT_EXPERIMENTDATA_H
#define PLOT_EXPERIMENTDATA_H

#include "plot_histogram/plot_histogram.h"

class plot_ExperimentData : public plot_histogram {
public:
    plot_ExperimentData();
    void run_inPixel();

private:
    static void set_2DSURFStyle(TCanvas* canvas, TH2D* hist);
    static TH2D* convert_toTH2D(TProfile2D* profile2D);

    std::vector<std::string> PIXEL_PITCH_, CHIP_TYPE_, VOLTAGE_, SEED_THRESHOLD_, NEIGHBOR_THRESHOLD_;
    std::string TIME_;
};

#endif // PLOT_EXPERIMENTDATA_H