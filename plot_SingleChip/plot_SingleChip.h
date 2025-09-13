#ifndef PLOT_SINGLECHIP_H
#define PLOT_SINGLECHIP_H

#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <cmath>
#include <algorithm>
#include <map>

// ROOT headers
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TColor.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TList.h>
#include <TObjString.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TCollection.h>

// My function headers
#include "Messenger/Messenger.h"
#include "plot_histogram/plot_histogram.h"
#include "plot_ExperimentData/plot_ExperimentData.h"

// Other function headers
#include "tools/json.hpp"
#include "tools/cxxopts.hpp"

enum class DataSource {
    onePixel,
    multiPixel
};

struct PlotHist1DType {
    std::string hist_name;
    std::string hist_path;
    std::string title;
    std::pair<double, double> x_range;
    std::pair<double, double> y_range;
    double scale_factor = 1.0;
};

struct PlotHist2DType {
    std::string hist_name;
    std::string hist_path;
    std::string title;
    std::pair<double, double> x_range;
    std::pair<double, double> y_range;
    std::pair<double, double> z_range;
    double scale_factor = 1.0;
};

struct PlotHist1DStyle {
    DataSource source;
    std::string pixel_pitch;
    std::string chip_type;
    std::string voltage;
    std::string carrior_model;
};

class plot_SingleChip {
public:
    plot_SingleChip();
    ~plot_SingleChip();

private:
    
    // class member parameters
    DataSource source_;
    
};

#endif