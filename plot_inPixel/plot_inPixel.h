#ifndef PLOT_INPIXEL_H
#define PLOT_INPIXEL_H

#include "plot_histogram/plot_histogram.h"

class plot_inPixel : public plot_histogram {
public:
    plot_inPixel();
    void run_inPixel();

private:
    std::vector<std::string> PIXEL_PITCH_, CHIP_TYPE_, VOLTAGE_, MODEL_, THRESHOLD_, PIXEL_NUMBER_;

    static TH2D* convert_toTH2D(TProfile2D* tprofile2D);
};

#endif // PLOT_INPIXEL_H