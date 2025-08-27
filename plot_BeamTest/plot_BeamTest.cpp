#include "plot_BeamTest.h"

plot_BeamTest::plot_BeamTest() {
    LOG_STATUS.source("plot_BeamTest::plot_BeamTest") << "plot_BeamTest object is created.";

    DATA_DIR_PATH_ = "/home/towa/alice3/hist/";
    TIME_ = plot_histogram::currentDateTime();
    canvas_ = new TCanvas("canvas", "canvas", 800, 600);
    gStyle->SetOptStat(0);

    title_latex_.SetTextSize(0.04);
    title_latex_.SetTextFont(62);
    condition_latex_.SetTextSize(0.03);
    condition_latex_.SetTextFont(62);
}

plot_BeamTest::~plot_BeamTest() {
}

void plot_BeamTest::run_plots(const std::vector<PlotConfig>& configs) {
    if(configs.empty()) {
        LOG_WARNING.source("plot_BeamTest::run_plots") << "No plot configuration provided.";
        return;
    }

    // resolution plot
    std::vector<TGraphErrors*> reso_graphs;
    for(const auto& conf : configs) {
        reso_graphs.push_back(create_graph_data(conf, "resolution"));
    }
    draw_multigraph("Resolution Comparison",
                    ";threshold [ADC];resolution in x [um]",
                    "plot/Combined_resolution.pdf",
                    configs,
                    reso_graphs,
                    //{2, 10.5},
                    {2, 13},
                    {40, 2510});
    //for(auto g : reso_graphs) delete g;

    // cluster size plot
    std::vector<TGraphErrors*> clsize_graphs;
    for(const auto& conf : configs) {
        clsize_graphs.push_back(create_graph_data(conf, "clustersize"));
    }
    draw_multigraph("Cluster Size Comparison",
                    ";threshold [ADC];mean cluster size",
                    "plot/Combined_ClusterSize.pdf",
                    configs,
                    clsize_graphs,
                    {1, 5.1},
                    {40, 2510});
    //for(auto g : clsize_graphs) delete g;

    // residual check
    LOG_STATUS.source("plot_BeamTest::run_plots") << "Creating residual check plots for all configurations.";
    for(const auto& config_to_check : configs) {
        std::string safe_label = config_to_check.legend_label;

        std::replace(safe_label.begin(), safe_label.end(), ' ', '_');
        std::replace(safe_label.begin(), safe_label.end(), '.', 'p');
        safe_label.erase(std::remove(safe_label.begin(), safe_label.end(), ','), safe_label.end());


        std::string output_filename = Form("plot/residuals/check_residuals_%s.pdf", safe_label.c_str());
        check_residual_fits(config_to_check, output_filename);
    }
}

void plot_BeamTest::draw_overlay_histograms(
    const std::string& canvas_title,
    const std::string& output_filename,
    const std::vector<PlotConfig>& configs,
    const std::string& neighbor_thd_for_all,
    const std::pair<double, double>& x_range
) {
    canvas_->Clear();
    canvas_->SetTopMargin(0.062);
    canvas_->SetBottomMargin(0.14);
    canvas_->SetLeftMargin(0.13);
    canvas_->SetRightMargin(0.07);

    TLegend* legend_point = new TLegend(0.50, 0.71, 0.90, 0.91);
    legend_point->SetFillStyle(0);
    legend_point->SetBorderSize(0);
    legend_point->SetTextSize(0.03);

    bool first_hist = true;
    std::vector<TH1D*> hists_to_draw;

    for(const auto& config : configs) {
        std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
                                           DATA_DIR_PATH_.c_str(),
                                           NAME_.c_str(),
                                           NAME_.c_str(),
                                           config.pixel_pitch.c_str(),
                                           config.chip_type.c_str(),
                                           config.voltage.c_str(),
                                           config.seed_thd.c_str(),
                                           neighbor_thd_for_all.c_str());
        
        TFile* inputFile = TFile::Open(input_file_path.c_str());
        if(!inputFile || inputFile->IsZombie()) {
            LOG_WARNING.source("plot_BeamTest::draw_overlay_histograms") << "Cannot open file " << input_file_path;
            if(inputFile) delete inputFile;
            continue;
        }

        TH1D* h_clcharge = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/cluster/clusterCharge", DUT_NAME_.c_str()));
        if(!h_clcharge) {
            LOG_WARNING.source("plot_BeamTest::draw_overlay_histograms") << "Histogram not found in " << input_file_path;
            inputFile->Close();
            delete inputFile;
            continue;
        }

        h_clcharge->SetDirectory(0);
        inputFile->Close();
        delete inputFile;

        h_clcharge->Rebin(10);
        if(h_clcharge->GetEntries() > 0) {
            h_clcharge->Scale(1.0 / h_clcharge->GetMaximum());
        }
        h_clcharge->SetLineColor(config.color);
        h_clcharge->SetMarkerColor(config.color);
        h_clcharge->SetMarkerStyle(config.marker_style);
        h_clcharge->SetLineWidth(2);

        legend_point->AddEntry(h_clcharge, config.legend_label.c_str(), "pe");

        if(first_hist) {
            h_clcharge->SetTitle(";charge [ADC];normalized counts");
            h_clcharge->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
            h_clcharge->GetYaxis()->SetRangeUser(0, 1.4);
            h_clcharge->Draw("PE");
        } else {
            h_clcharge->Draw("samePE");
        }

        legend_point->Draw();
        condition_latex_.DrawLatexNDC(0.15, 0.91, BEAM_INFO_.c_str());
        condition_latex_.DrawLatexNDC(0.15, 0.87, Form("Neighbor Threshold = %s ADC", neighbor_thd_for_all.c_str()));
        condition_latex_.DrawLatexNDC(0.15, 0.83, Form("Plotted on %s", TIME_.c_str()));

        canvas_->SaveAs(output_filename.c_str());
    }
}

TGraphErrors* plot_BeamTest::create_graph_data(
    const PlotConfig& config,
    const std::string& quantity_to_extract) {
    std::string name, dut_name;
    if(config.source == DataSource::KEK202412) {
        NAME_ = "kek202412";
        DUT_NAME_ = "CE65_3";
    } else if(config.source == DataSource::SPS202404) {
        NAME_ = "sps202404";
        DUT_NAME_ = "CE65_6";
    }
    
    std::vector<double> x_vals, y_vals, x_errs, y_errs;

    for(const auto& thd : config.scan_values) {
        double neighbor_val = std::stod(thd);
        double seed_val = std::stod(config.seed_thd);

        std::string seed_thd_for_file = config.seed_thd;
        const std::string& neighbor_thd_for_file = thd;
        if(neighbor_val > seed_val) {
            seed_thd_for_file = thd;
        }

        std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
                                           DATA_DIR_PATH_.c_str(),
                                           NAME_.c_str(),
                                           NAME_.c_str(),
                                           config.pixel_pitch.c_str(),
                                           config.chip_type.c_str(),
                                           config.voltage.c_str(),
                                           seed_thd_for_file.c_str(), //config.seed_thd.c_str(),
                                           thd.c_str());

        TFile* inputFile = TFile::Open(input_file_path.c_str());
        if(!inputFile || inputFile->IsZombie()) {
            LOG_WARNING.source("plot_BeamTest::created_graph_data") << "Cannot open file " << input_file_path;
            if(inputFile) delete inputFile;
                continue;
        }

        TH1D* h_residual = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str()));
        TH1D* h_clsize = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/cluster/clusterSize", DUT_NAME_.c_str()));

        if(!h_residual || !h_clsize) {
            LOG_WARNING.source("plot_BeamTest::create_graph_data") << "Histogram not found in " << input_file_path;
            inputFile->Close();
            delete inputFile;
            continue;
        }

        double y_val = 0.0;
        double y_err = 0.0;

        if(quantity_to_extract == "resolution") {
            h_residual->Rebin(2);
            TF1* fResidual = plot_histogram::optimise_hist_gaus(h_residual, kBlack);
            y_val = fResidual->GetParameter(2);
            y_err = fResidual->GetParError(2);
        } else if (quantity_to_extract == "clustersize") {
            y_val = h_clsize->GetMean();
            y_err = h_clsize->GetMeanError();
        }

        x_vals.push_back(std::stod(thd));
        y_vals.push_back(y_val);
        x_errs.push_back(0);
        y_errs.push_back(y_err);

        inputFile->Close();
        delete inputFile;
    }

    return new TGraphErrors(x_vals.size(), x_vals.data(), y_vals.data(), x_errs.data(), y_errs.data());
}

void plot_BeamTest::draw_multigraph(
    const std::string& canvas_title,
    const std::string& title,
    const std::string& output_filename,
    const std::vector<PlotConfig>& configs,
    const std::vector<TGraphErrors*>& graphs,
    const std::pair<double, double>& y_range,
    const std::pair<double, double>& x_range) {
    canvas_->Clear();
    canvas_->SetTopMargin(0.062);
    canvas_->SetBottomMargin(0.14);
    canvas_->SetLeftMargin(0.13);
    canvas_->SetRightMargin(0.07);
    canvas_->SetGrid();
    gStyle->SetGridStyle(1);
    gStyle->SetGridColor(kGray);

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(title.c_str());

    double y_size = graphs.size() * 0.02;
    double legend_y_min = 0.91 - y_size;

    TLegend* legend_point = new TLegend(0.40, legend_y_min, 0.90, 0.91);
    legend_point->SetFillStyle(0);
    legend_point->SetBorderSize(0);
    legend_point->SetTextSize(0.03);
    legend_point->SetNColumns(2);

    TLegend* legend_line = new TLegend(0.40, legend_y_min, 0.90, 0.91);
    legend_line->SetFillStyle(0);
    legend_line->SetBorderSize(0);
    legend_line->SetTextSize(0.03);
    legend_line->SetTextColor(kWhite);
    legend_line->SetNColumns(2);

    for(size_t i=0; i<configs.size(); ++i) {
        if(!graphs[i] || graphs[i]->GetN() == 0 || graphs[i] == nullptr) {
            LOG_WARNING.source("plot_BeamTest::draw_multigraph") << "Graph not found.";
            continue; 
        }

        //graphs[i]->SetLineColorAlpha(configs[i].color, 0.5);
        graphs[i]->SetMarkerColor(configs[i].color);
        graphs[i]->SetMarkerStyle(configs[i].marker_style);
        graphs[i]->SetLineStyle(configs[i].line_style);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerSize(configs[i].marker_size);

        TGraphErrors* graph_point = (TGraphErrors*)graphs[i]->Clone("graph_point");
        graph_point->SetLineStyle(1);

        graphs[i]->SetLineColorAlpha(configs[i].color, 0.5);
        graph_point->SetLineColor(configs[i].color);

        mg->Add(graphs[i], "L");
        mg->Add(graph_point, "PE");
        legend_line->AddEntry(graphs[i], configs[i].legend_label.c_str(), "l");
        legend_point->AddEntry(graph_point, configs[i].legend_label.c_str(), "pe");
    }

    if(mg->GetListOfGraphs() == nullptr) {
        LOG_ERROR.source("plot_BeamTest::draw_multigraph") << "No valid graphs to draw for " << canvas_title;
        delete mg;
        delete legend_line;
        delete legend_point;
        return;
    }

    mg->Draw("A");
    mg->GetYaxis()->SetRangeUser(y_range.first, y_range.second);
    mg->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetLabelSize(0.04);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetTitleOffset(0.8);

    legend_line->Draw();
    legend_point->Draw();
    //condition_latex_.DrawLatexNDC(0.15, 0.91, BEAM_INFO_.c_str());
    condition_latex_.DrawLatexNDC(0.68, 0.91, Form("Plotted on %s", TIME_.c_str()));
    //title_latex_.DrawLatexNDC(0.15, 0.83, canvas_title.c_str());

    canvas_->SaveAs(output_filename.c_str());
}

void plot_BeamTest::check_residual_fits(
    const PlotConfig& config,
    const std::string& output_filename) {
    LOG_STATUS.source("plot_BeamTest::check_residual_fits") << "Creating residual check plots for " << config.legend_label;

    //std::string name, dut_name, beam_info;
    if (config.source == DataSource::KEK202412) {
        NAME_ = "kek202412";
        DUT_NAME_ = "CE65_3";
        BEAM_INFO_ = "e^{-} 3GeV/c @KEK-PFAR (Dec. 2024)";
    } else { // DataSource::SPS202404
        NAME_ = "sps202404";
        DUT_NAME_ = "CE65_6";
        BEAM_INFO_ = "hadron 120GeV/c @CERN-SPS (Apr. 2024)";
    }
    const std::string DATA_DIR_PATH_ = "/home/towa/alice3/hist/";


    const int max_pads_per_canvas = 12;
    const int n_cols = 4;
    const int n_rows = 3;
    const std::vector<std::string>& thresholds = config.scan_values;
    const int n_total_pads = thresholds.size();
    
    if(n_total_pads == 0) {
        LOG_WARNING.source("plot_BeamTest::check_residual_fits") << "No scan values provided for residual check.";
        return;
    }

    TCanvas* check_canvas = new TCanvas("check_canvas", "residual fits", 1920, 1080);
    check_canvas->Print(Form("%s[", output_filename.c_str())); // Open pdf

    // title page
    check_canvas->Clear();
    TLatex title_page_latex;
    title_page_latex.SetTextAlign(12);

    title_page_latex.SetTextSize(0.05);
    title_page_latex.DrawLatexNDC(0.1, 0.85, "Residual Distribution & Fitting Results");

    TLatex information_latex;
    information_latex.SetTextSize(0.04);
    information_latex.DrawLatexNDC(0.1, 0.75, Form("Condition: %s", config.legend_label.c_str()));
    information_latex.DrawLatexNDC(0.1, 0.70, Form("Experiment: %s", BEAM_INFO_.c_str()));
    information_latex.DrawLatexNDC(0.1, 0.65, Form("Seed Threshold: %s ADC", config.seed_thd.c_str()));
    
    information_latex.DrawLatexNDC(0.1, 0.55, "Neighbor Thresholds Scanned [ADC]:");
    
    information_latex.SetTextSize(0.035);
    double y_pos = 0.50;
    std::string current_line;
    for(size_t i = 0; i < thresholds.size(); ++i) {
        current_line += thresholds[i];
        if (i < thresholds.size() - 1) {
            current_line += ", ";
        }

        if (current_line.length() > 60 || i == thresholds.size() - 1) {
            title_page_latex.DrawLatexNDC(0.15, y_pos, current_line.c_str());
            y_pos -= 0.04;
            current_line.clear();
        }
    }

    check_canvas->Print(output_filename.c_str());

    // plot page
    int pads_drawn = 0;
    while(pads_drawn < n_total_pads) {
        check_canvas->Clear();

        int pads_on_this_canvas = std::min(max_pads_per_canvas, n_total_pads - pads_drawn);

        // int n_cols = (pads_on_this_canvas > 9) ? 4 : 3;
        // int n_rows = static_cast<int>(std::ceil(static_cast<double>(pads_on_this_canvas) / n_cols));
        check_canvas->Divide(n_cols, n_rows);

        for(int i=0; i<pads_on_this_canvas; ++i) {
            int current_index = pads_drawn + i;
            check_canvas->cd(i + 1)->SetMargin(0.15, 0.05, 0.12, 0.05);

            const std::string& thd = thresholds[current_index];
            double neighbor_val = std::stod(thd);
            double seed_val = std::stod(config.seed_thd);

            std::string seed_thd_for_file = config.seed_thd;
            //const std::string7 neighbor_thd_for_file = thd;
            if(neighbor_val >  seed_val) {
                seed_thd_for_file = thd;
            }

            //const std::string& thd = thresholds[current_index];
            std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
                                               DATA_DIR_PATH_.c_str(),
                                               NAME_.c_str(),
                                               NAME_.c_str(),
                                               config.pixel_pitch.c_str(),
                                               config.chip_type.c_str(),
                                               config.voltage.c_str(),
                                               seed_thd_for_file.c_str(),  //config.seed_thd.c_str(),
                                               thd.c_str());
            
            TFile* inputFile = TFile::Open(input_file_path.c_str());
            if(!inputFile || inputFile->IsZombie()) {
                LOG_WARNING.source("plot_BeamTest::check_residual_fits") << "Cannot open file " << input_file_path;
                if(inputFile) delete inputFile;
                continue;
            }

            TH1D* h_res = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str()));
            if(!h_res) {
                LOG_WARNING.source("plot_BeamTest::check_residual_fits") << "Histogram not found in " << input_file_path;
                delete inputFile;
                continue;
            }

            h_res->SetDirectory(0);
            inputFile->Close();
            delete inputFile;

            h_res->SetTitle(";x_{track}-x_{hit} [um];normalized counts");
            //h_res->GetXaxis()->SetRangeUser(-50, 50);
            h_res->GetXaxis()->SetTitleSize(0.06);
            h_res->GetYaxis()->SetTitleSize(0.06);
            h_res->GetXaxis()->SetLabelSize(0.05);
            h_res->GetYaxis()->SetLabelSize(0.05);
            h_res->Rebin(2);

            h_res->Draw("PE");

            // fitting
            // TF1* fit_res = new TF1("fit_res", "gaus", -50, 50);
            // fit_res->SetLineColor(kRed);
            // fit_res->SetLineWidth(2);
            // h_res->Fit(fit_res, "NSLQ+", "", -20, 20);
            TF1* fit_res = plot_histogram::optimise_hist_gaus(h_res, kRed);
            fit_res->SetLineWidth(2);
            h_res->GetXaxis()->SetRangeUser(-50, 50);

            TLatex latex_parameter;
            latex_parameter.SetTextSize(0.05);
            latex_parameter.SetTextColor(kBlack);
            latex_parameter.DrawLatexNDC(0.17, 0.9, Form("#sigma = %.2f #pm %.2f um", fit_res->GetParameter(2), fit_res->GetParError(2))); 
            latex_parameter.DrawLatexNDC(0.17, 0.84, Form("NeighborThd = %s ADC", thd.c_str()));          
        }

        check_canvas->Print(output_filename.c_str());
        pads_drawn += pads_on_this_canvas;
    }

    check_canvas->Print(Form("%s]", output_filename.c_str()));
    delete check_canvas;
}

Color_t plot_BeamTest::string_to_ROOTColor(const std::string& color_str) {
    static const std::map<std::string, Color_t> color_map = {
        {"kBlack", kBlack},
        {"kWhite", kWhite},
        {"kGray", kGray},
        {"kRed", kRed},
        {"kRed+1", kRed + 1},
        {"kRed+2", kRed + 2},
        {"kRed+3", kRed + 3},
        {"kGreen", kGreen},
        {"kGreen+1", kGreen + 1},
        {"kGreen+2", kGreen + 2},
        {"kBlue", kBlue},
        {"kBlue+1", kBlue + 1},
        {"kBlue+2", kBlue + 2},
        {"kYellow", kYellow},
        {"kMagenta", kMagenta},
        {"kCyan", kCyan},
        {"kAzure+2", kAzure + 2},
        {"kPink+7", kPink + 7}
    };

    auto it = color_map.find(color_str);
    if(it != color_map.end()) {
        return it->second;
    }

    LOG_WARNING.source("plot_BeamTest::string_to_ROOTColor") << "Color string '" << color_str << "' not found. Using kBlack.";
    return kBlack;
}

std::vector<PlotConfig> plot_BeamTest::load_jsonConfigs(const std::string& filename) {
    std::vector<PlotConfig> configs;
    std::ifstream file(filename);
    if(!file.is_open()) {
        LOG_ERROR.source("plot_BeamTest::load_jsonConfigs") << "Could not open config file: " << filename;
        return configs;
    }

    nlohmann::json j;
    try {
        file >> j;
    } catch (nlohmann::json::parse_error& e) {
        LOG_ERROR.source("plot_BeamTest::load_jsonConfigs") << "JSON parse error: " << e.what();
        return configs;
    }

    for(const auto& item : j) {
        PlotConfig conf;

        // Conversion to DataSource enum
        std::string source_str = item.at("source");
        if(source_str == "KEK202412") {
            conf.source = DataSource::KEK202412;
        } else if (source_str == "SPS202404") {
            conf.source = DataSource::SPS202404;
        } else {
            LOG_WARNING.source("plot_BeamTest::load_jsonConfigs") << "Unknown data source: " << source_str;
            continue;
        }

        conf.pixel_pitch    = item.at("pixel_pitch").get<std::string>();
        conf.chip_type      = item.at("chip_type").get<std::string>();
        conf.voltage        = item.at("voltage").get<std::string>();
        conf.seed_thd       = item.at("seed_thd").get<std::string>();
        conf.scan_values    = item.at("scan_values").get<std::vector<std::string>>();
        conf.legend_label   = item.at("legend_label").get<std::string>();
        //conf.color          = item.at("color").get<std::string>();
        conf.marker_style   = item.at("marker_style").get<int>();
        conf.marker_size    = item.at("marker_size").get<double>();
        conf.line_style     = item.at("line_style").get<int>();

        std::string color_string  = item.at("color").get<std::string>();
        conf.color          = string_to_ROOTColor(color_string);

        configs.push_back(conf);
    }

    return configs;
}

void plot_BeamTest::BeamTest_main(int argc, char* argv[]) {
    cxxopts::Options options("plot_BeamTest::BeamTest_main", "Beamtest object code");
    options.add_options()
        ("f,file", "Json file name", cxxopts::value<std::string>())
        ("h,help", "show help message");
    
    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        std::cout << options.help() << std::endl;
        return;
    }

    std::string config_filename;
    if(result.count("file")) {
        config_filename = result["file"].as<std::string>();
    } else {
        LOG_WARNING.source("plot_BeamTest::BeamTest_main") << "JSON file name is NOT defined. Using default file name: plot_BeamTest/json/analysis.json";
        config_filename = "plot_BeamTest/json/analysis.json";
    }

    std::vector<PlotConfig> my_comparison = plot_BeamTest::load_jsonConfigs(config_filename);

    if(my_comparison.empty()) {
        LOG_ERROR.source("plot_BeamTest::BeamTest_main") << "No configurations were loaded.";
        return;
    }

    plot_BeamTest::run_plots(my_comparison);
}