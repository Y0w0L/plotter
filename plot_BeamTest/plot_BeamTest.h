#ifndef PLOT_BEAMTEST_H
#define PLOT_BEAMTEST_H

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

// My function headers
#include "Messenger/Messenger.h"
#include "plot_histogram/plot_histogram.h"

// Other function headers
#include "tools/json.hpp"
#include "tools/cxxopts.hpp"

// Data source
enum class DataSource {
    KEK202412,
    SPS202404
};

struct PlotConfig {
    DataSource source;
    std::string pixel_pitch;
    std::string chip_type;
    std::string voltage;
    std::string seed_thd;
    std::vector<std::string> scan_values;

    // style information
    std::string legend_label;
    int color = kBlack;
    int marker_style = 20;
    double marker_size = 1.2;
    int line_style = 1;
};

class plot_BeamTest {
public:
    // Assign data source in constructer
    //plot_BeamTest(DataSource source);
    plot_BeamTest();
    ~plot_BeamTest();

    // plot for each data surce
    void run_kek_plots();
    void run_sps_plots();

    void run_plots(const std::vector<PlotConfig>& configs);

    void check_residual_fits(
        const PlotConfig& config,
        const std::string& quantity_to_extract
    );

    static std::vector<PlotConfig> load_jsonConfigs(const std::string& filename);

    void BeamTest_main(int argc, char* argv[]);

private:
    // helper functions
    TGraphErrors* create_graph_data(
        const PlotConfig& config,
        const std::string& quantity_to_extract // "resolution" or "clustersize"
    );

    void draw_multigraph(
        const std::string& canvas_title,
        const std::string& title, // 後で変えるかも
        const std::string& output_filename,
        const std::vector<PlotConfig>& configs,
        const std::vector<TGraphErrors*>& graphs,
        const std::pair<double, double>& y_range,
        const std::pair<double, double>& x_range
    );

    void draw_overlay_histograms(
        const std::string& canvas_title,
        const std::string& output_filename,
        const std::vector<PlotConfig>& configs,
        const std::string& neighbor_thd_for_all,
        const std::pair<double, double>& x_range
    );

    static Color_t string_to_ROOTColor(const std::string& color_str);

    // void check_residual_fits(
    //     const PlotConfig& config,
    //     const std::string& output_filename
    // );

    // class member parameters
    DataSource source_;
    std::string TIME_;
    std::string NAME_;
    std::string DUT_NAME_;
    std::string DATA_DIR_PATH_;
    std::string BEAM_INFO_;
    TCanvas* canvas_;
    TLatex title_latex_;
    TLatex condition_latex_;
};

#endif // PLOT_BEAMTEST_H