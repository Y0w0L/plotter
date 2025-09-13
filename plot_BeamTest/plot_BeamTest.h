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

// Data source
enum class DataSource {
    KEK202412,
    SPS202404,
    SingleChipSim,
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

struct ChipParameters {
    std::string pixel_pitch;
    std::string chip_type;
    std::string voltage;
    std::string seed_thd;
    std::string neighbor_thd;
};

struct InPixelPlotConfig {
    std::string name;
    std::string hist_path;
    std::string title;
    std::string z_axis_title;
    double z_min;
    double z_max;
    double scale_factor = 1.0;
};

struct PathPoint {
    std::string label;
    std::string x_expr; // x coordinates formula (ex, "half_pitch - inset_x")
    std::string y_expr;
};

struct PathSegment {
    std::string from; // Start point
    std::string to;   // Goal point
    std::string label;
    std::string color_str;
};

struct PathConfig {
    std::string name;
    std::vector<PathPoint> points;
    std::vector<PathSegment> segments;
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

    void BeamTest_main(int argc, char* argv[]);

    void run_inPixelAnalysis(const std::vector<PlotConfig>& configs, const std::vector<InPixelPlotConfig>& plot_types);

    void run_inPixelPathAnalysis(const std::vector<PlotConfig>& configs, const std::vector<InPixelPlotConfig>& plot_types);

    // void run_inPixelPathAnalysis(
    //     const std::vector<PlotConfig>& plotConfigs,
    //     const std::vector<InPixelPlotConfig>& inPixelPlotConfigs,
    //     const std::vector<PathConfig>& pathConfigs
    // );

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

    void drawBeamInfo(const std::string& beam_info, double x, double y);
    void drawChipInfo(const ChipParameters& params, double x, double y);
    void createAndSaveInPixelPlot(
        TH2D* hist,
        const InPixelPlotConfig& config,
        const ChipParameters& params,
        const std::string& beam_info,
        TDirectory* output_dir,
        const std::string& base_path,
        const std::string& chip_variation_name
    );

    double evaluate_expr(const std::string&  expr, double half_pitch, double inset_x, double inset_y);

    std::vector<PlotConfig> load_jsonConfigs(const nlohmann::json& j);
    std::vector<InPixelPlotConfig> load_jsonInPixelPlotConfigs(const nlohmann::json& j);
    std::vector<PathConfig> load_jsonPathConfigs(const nlohmann::json& j);

    template<typename THist>
    THist* get_merged_object(const std::string& base_file_path, const std::string& object_name);

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