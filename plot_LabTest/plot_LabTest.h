#ifndef PLOT_LABTEST_H
#define PLOT_LABTEST_H

#include "plot_histogram/plot_histogram.h"

class plot_LabTest : public plot_histogram {
public:
    plot_LabTest();

    static void set_TH1DStyle(TH1D* hist);
    static void set_TH1DPosition(TCanvas* canvas);
    void run_LabTest();

private:
    std::string FILENAME_STD_1V, FILENAME_STD_4V, FILENAME_STD_7V, FILENAME_STD_10V;
    std::string FILENAME_BLK_1V, FILENAME_BLK_4V, FILENAME_BLK_7V, FILENAME_BLK_10V;
    std::string FILENAME_GAP_1V, FILENAME_GAP_4V, FILENAME_GAP_7V, FILENAME_GAP_10V;
    std::vector<std::string> FILENAMES_STD_, FILENAMES_BLK_, FILENAMES_GAP_;
    std::vector<std::vector<std::string>> FILENAMES_;
    std::vector<std::string> CHIP_TYPE_, VOLTAGE_;

    std::string TIME_;
};

#endif // PLOT_LABTEST_H