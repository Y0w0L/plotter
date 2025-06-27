#ifndef PLOT_INPIXEL_H
#define PLOT_INPIXEL_H

#include "plot_histogram/plot_histogram.h"

class plot_inPixel : public plot_histogram {
public:
    plot_inPixel();
    void run_inPixel();

private:
    std::vector<std::string> PIXEL_PITCH_, CHIP_TYPE_, VOLTAGE_, MODEL_, THRESHOLD_, PIXEL_NUMBER_;
    std::string TIME_;

    static void set_2DSURFStyle(TCanvas* canvas, TH2D* hist);
    static TH2D* convert_toTH2D(TProfile2D* tprofile2D);
    static std::string replace_underscore_to_slash(std::string text);
};

#endif // PLOT_INPIXEL_H