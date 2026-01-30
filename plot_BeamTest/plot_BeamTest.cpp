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

    use_electron_scale_ = false;
    extract_tracking_resolution_ = false;
}

plot_BeamTest::~plot_BeamTest() {
}

template<typename THist>
THist* plot_BeamTest::get_merged_object(
    const std::string& base_file_path,
    const std::string& object_name
) {
    TString dir_path = gSystem->DirName(base_file_path.c_str());
    TString base_name = gSystem->BaseName(base_file_path.c_str());

    void* drip = gSystem->OpenDirectory(dir_path);
    if(!drip) {
        LOG_WARNING.source("plot_BeamTest::get_merged_object") << "Cannot open directory: " << dir_path.Data();
        return nullptr;
    }

    std::vector<std::string> matching_files;
    const char* entry;
    while((entry = gSystem->GetDirEntry(drip))) {
        std::string filename = entry;
        // search for filename=base_name + ~ + .root
        if(filename.rfind(base_name.Data(), 0) == 0 && filename.find(".root") != std::string::npos) {
            matching_files.push_back(std::string(dir_path.Data()) + "/" + filename);
        }
    }
    gSystem->FreeDirectory(drip);

    if(matching_files.empty()) {
        LOG_WARNING.source("plot_BeamTest::get_merged_object") << "No matching files found for base: " << base_file_path;
        return nullptr;
    }

    LOG_DEBUG.source("plot_BeamTest::get_merged_object") << "Found " << matching_files.size() << " files for merging for base: " << base_name.Data();
    
    THist* merged_obj = nullptr;
    for(const auto& file_path : matching_files) {
        //std::cout << "----------------" << std::endl;
        TFile* file = TFile::Open(file_path.c_str());
        if(!file || file->IsZombie()) {
            LOG_WARNING.source("plot_BeamTest::get_merged_object") << "Cannot open file " << file_path;
            if(file) delete file;
            continue;
        }

        THist* h = (THist*)file->Get(object_name.c_str());
        if(!h) {
            LOG_WARNING.source("plot_BeamTest::get_merged_object") << "Object '" << object_name << "' not found in " << file_path;
            file->Close();
            delete file;
            continue;
        }

        if(!merged_obj) {
            merged_obj = (THist*)h->Clone();
            merged_obj->SetDirectory(0);
        } else {
            //std::cout << "====================" << std::endl;
            merged_obj->Add(h);
        }

        file->Close();
        delete file;
    }

    return merged_obj;
}

void plot_BeamTest::run_plots(const std::vector<PlotConfig>& configs) {
    if(configs.empty()) {
        LOG_WARNING.source("plot_BeamTest::run_plots") << "No plot configuration provided.";
        return;
    }

    std::string x_axis_label = use_electron_scale_ ? ";threshold [e^{-}]" : ";threshold [ADC]";

    // resolution plot
    std::vector<TGraphErrors*> reso_graphs;
    for(const auto& conf : configs) {
        reso_graphs.push_back(create_graph_data(conf, "resolution"));
    }

    std::optional<std::pair<double, double>> reso_x_range = std::nullopt;
    if(!use_electron_scale_) {
        reso_x_range = {0, 2510};
    }

    std::string target_thd = "200";

    // Cluster Charge の重ね書き
    draw_overlay_histograms(
        "Cluster Charge Comparison", 
        Form("plot/Overlay_ClusterCharge_Thd%s.pdf", target_thd.c_str()), 
        configs, target_thd, {0, 6000}, "clusterCharge"
    );

    // Seed Charge の重ね書き
    draw_overlay_histograms(
        "Seed Charge Comparison", 
        Form("plot/Overlay_SeedCharge_Thd%s.pdf", target_thd.c_str()), 
        configs, target_thd, {0, 4000}, "seedCharge"
    );

    draw_multigraph("Resolution Comparison",
                    (x_axis_label + ";RMS [um]").c_str(),
                    "plot/Combined_resolution.pdf",
                    configs,
                    reso_graphs,
                    //{2, 10.5},
                    {2.7, 8.1},
                    reso_x_range);
    //for(auto g : reso_graphs) delete g;

    // cluster size plot
    std::vector<TGraphErrors*> clsize_graphs;
    for(const auto& conf : configs) {
        clsize_graphs.push_back(create_graph_data(conf, "clustersize"));
    }

    std::optional<std::pair<double, double>> clsize_x_range = std::nullopt;
    if(!use_electron_scale_) {
        clsize_x_range = {0, 2510};
    }

    draw_multigraph("Cluster Size Comparison",
                    (x_axis_label + ";mean cluster size").c_str(),
                    "plot/Combined_ClusterSize.pdf",
                    configs,
                    clsize_graphs,
                    {0.9, 4.55},
                    clsize_x_range);
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
    const std::pair<double, double>& x_range,
    const std::string& target_quantity
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
        if(config.source == DataSource::KEK202412) {
            NAME_ = "kek202412";
            DUT_NAME_ = "CE65_3";
            BEAM_INFO_ = "@KEK PF-AR Dec. 2024, 3 GeV/c electrons";
        }else if(config.source == DataSource::SPS202404) {
            NAME_ = "sps202404";
            DUT_NAME_ = "CE65_6";
            BEAM_INFO_ = "@CERN SPS Apr. 2024, 120 GeV/c hadrons";
        }else if(config.source == DataSource::SingleChipSim) {
            NAME_ = "ce65sim202505";
            DUT_NAME_ = "CE65";
            BEAM_INFO_ = "120 GeV/c pions";
        } else if(config.source == DataSource::SingleChipDrift) {
            NAME_ = "ce65driftTime";
            DUT_NAME_ = "CE65";
            BEAM_INFO_ = "120 GeV/c pions";
        } else {
            LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "No Datasource!";
            return;
        }

        std::string input_file_path;
        double seed_val = std::stod(config.seed_thd);
        double neighbor_val = std::stod(neighbor_thd_for_all);
        std::string seed_thd_for_file = (neighbor_val > seed_val) ? std::to_string(neighbor_val) : config.seed_thd;

        if (!config.base_file_name.empty()) {
            // JSONで指定がある場合
            if (config.base_file_name.find("%s") != std::string::npos) {
                // neighbor_thd_for_all を引数としてフォーマット
                input_file_path = Form(config.base_file_name.c_str(), seed_thd_for_file.c_str(), neighbor_thd_for_all.c_str());
            } else {
                input_file_path = config.base_file_name;
            }
        } else {
            // 既存ロジック
            if(config.source == DataSource::SingleChipSim || config.source == DataSource::SingleChipDrift) {
                seed_thd_for_file = "0";
            }

            input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
                                   DATA_DIR_PATH_.c_str(),
                                   NAME_.c_str(),
                                   NAME_.c_str(),
                                   config.pixel_pitch.c_str(),
                                   config.chip_type.c_str(),
                                   config.voltage.c_str(),
                                   seed_thd_for_file.c_str(),
                                   neighbor_thd_for_all.c_str());
        }

        // double seed_val = std::stod(config.seed_thd);
        // double neighbor_val = std::stod(neighbor_thd_for_all);
        // std::string seed_thd_for_file = (neighbor_val > seed_val) ? std::to_string(neighbor_val) : config.seed_thd;
        // if(config.source == DataSource::SingleChipSim || config.source == DataSource::SingleChipDrift) {
        //     seed_thd_for_file = "0";
        // }

        // std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
        //                                    DATA_DIR_PATH_.c_str(),
        //                                    NAME_.c_str(),
        //                                    NAME_.c_str(),
        //                                    config.pixel_pitch.c_str(),
        //                                    config.chip_type.c_str(),
        //                                    config.voltage.c_str(),
        //                                    seed_thd_for_file.c_str(),
        //                                    neighbor_thd_for_all.c_str());

        // std::string hist_name = Form("AnalysisCE65/%s/cluster/clusterCharge", DUT_NAME_.c_str());
        // if(config.source == DataSource::SingleChipSim) {
        //     hist_name = "cluster_charge";
        // }

        // std::string hist_name;
        // if (!config.hist_path.empty()) {
        //     hist_name = config.hist_path;
        // } else {
        //     // 既存ロジック
        //     hist_name = Form("AnalysisCE65/%s/cluster/clusterCharge", DUT_NAME_.c_str());
        //     if(config.source == DataSource::SingleChipSim) {
        //         hist_name = "cluster_charge";
        //     }
        // }
        std::string hist_name;
        if (target_quantity == "seedCharge") {
            hist_name = Form("AnalysisCE65/%s/cluster/clusterSeedCharge", DUT_NAME_.c_str());
            if(config.source == DataSource::SingleChipSim) hist_name = "seed_charge";
        } else {
            // デフォルトは clusterCharge
            hist_name = !config.hist_path.empty() ? config.hist_path : Form("AnalysisCE65/%s/cluster/clusterCharge", DUT_NAME_.c_str());
            if(config.source == DataSource::SingleChipSim) hist_name = "cluster_charge";
        }

        TH1D* h_clcharge = get_merged_object<TH1D>(input_file_path, hist_name);
        if(!h_clcharge) {
            LOG_WARNING.source("plot_BeamTest::draw_overlay_histograms") << "Merged histogram could not be created for base " << input_file_path;
            continue;
        }

        h_clcharge->SetDirectory(0);
        // inputFile->Close();
        // delete inputFile;

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
    double tracking_resolution;
    if(config.source == DataSource::KEK202412) {
        NAME_ = "kek202412";
        DUT_NAME_ = "CE65_3";
        tracking_resolution = 0;
    } else if(config.source == DataSource::SPS202404) {
        NAME_ = "sps202404";
        DUT_NAME_ = "CE65_6";
        tracking_resolution = 2.85;
    } else if(config.source == DataSource::SingleChipSim) {
        NAME_ = "ce65sim202505";
        DUT_NAME_ = "CE65";
        tracking_resolution = 0;
    } else {
        LOG_ERROR.source("plot_BeamTest::create_graph_data");
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

        // if(config.source == DataSource::SingleChipSim) {
        //     seed_thd_for_file = "0";
        // }

        // std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
        //                                    DATA_DIR_PATH_.c_str(),
        //                                    NAME_.c_str(),
        //                                    NAME_.c_str(),
        //                                    config.pixel_pitch.c_str(),
        //                                    config.chip_type.c_str(),
        //                                    config.voltage.c_str(),
        //                                    seed_thd_for_file.c_str(), //config.seed_thd.c_str(),
        //                                    thd.c_str());

        // TFile* inputFile = TFile::Open(input_file_path.c_str());
        // if(!inputFile || inputFile->IsZombie()) {
        //     LOG_WARNING.source("plot_BeamTest::created_graph_data") << "Cannot open file " << input_file_path;
        //     if(inputFile) delete inputFile;
        //         continue;
        // }

        // TH1D* h_residual = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str()));
        // TH1D* h_clsize = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/cluster/clusterSize", DUT_NAME_.c_str()));

        // if(!h_residual || !h_clsize) {
        //     LOG_WARNING.source("plot_BeamTest::create_graph_data") << "Histogram not found in " << input_file_path;
        //     inputFile->Close();
        //     delete inputFile;
        //     continue;
        // }

        // std::string base_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
        //                                    DATA_DIR_PATH_.c_str(),
        //                                    NAME_.c_str(),
        //                                    NAME_.c_str(),
        //                                    config.pixel_pitch.c_str(),
        //                                    config.chip_type.c_str(),
        //                                    config.voltage.c_str(),
        //                                    seed_thd_for_file.c_str(),
        //                                    thd.c_str());

        // std::string hist_name_residual = Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str());
        // std::string hist_name_clsize = Form("AnalysisCE65/%s/cluster/clusterSize", DUT_NAME_.c_str());

        // if(config.source == DataSource::SingleChipSim) {
        //     hist_name_residual = "residual_x";
        //     hist_name_clsize = "cluster_size";
        // }

        std::string base_file_path;

        if(!config.base_file_name.empty()) {
            if(config.base_file_name.find("%s") != std::string::npos) {
                base_file_path = Form(config.base_file_name.c_str(), seed_thd_for_file.c_str(), thd.c_str());
            } else {
                base_file_path = config.base_file_name;
            }
        } else {
            if(config.source == DataSource::SingleChipSim || config.source == DataSource::SingleChipDrift) {
                seed_thd_for_file = "0";
            }
            
            base_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
                                   DATA_DIR_PATH_.c_str(),
                                   NAME_.c_str(),
                                   NAME_.c_str(),
                                   config.pixel_pitch.c_str(),
                                   config.chip_type.c_str(),
                                   config.voltage.c_str(),
                                   seed_thd_for_file.c_str(),
                                   thd.c_str());
        }

        std::string hist_name_residual;
        std::string hist_name_clsize;

        if(!config.hist_path_residual.empty()) {
            hist_name_residual = config.hist_path_residual;
        } else {
            hist_name_residual = Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str());
            if(config.source == DataSource::SingleChipSim) hist_name_residual = "residual_x";
        }
        if (!config.hist_path_clsize.empty()) {
            hist_name_clsize = config.hist_path_clsize;
        } else {
            // デフォルト
            hist_name_clsize = Form("AnalysisCE65/%s/cluster/clusterSize", DUT_NAME_.c_str());
            if(config.source == DataSource::SingleChipSim) hist_name_clsize = "cluster_size";
        }

        TH1D* h_residual = get_merged_object<TH1D>(base_file_path, hist_name_residual);
        TH1D* h_clsize = get_merged_object<TH1D>(base_file_path, hist_name_clsize);

        if(!h_residual || !h_clsize) {
            LOG_WARNING.source("plot_BeamTest::create_graph_data") << "One or more merged histograms not found for base " << base_file_path;
            if(h_residual) delete h_residual;
            if(h_clsize) delete h_clsize;
            continue;
        }

        double y_val = 0.0;
        double y_err = 0.0;
        double y_rms = 0.0;

        double x_val = std::stod(thd);
        bool is_exp_data = (config.source == DataSource::KEK202412 || config.source == DataSource::SPS202404);

        if(use_electron_scale_ && is_exp_data && config.adc_to_electron_factor > 0.0) {
            x_val *= config.adc_to_electron_factor;
        }

        if(quantity_to_extract == "resolution") {
            h_residual->Rebin(4);
            TF1* fResidual = plot_histogram::optimise_hist_gaus(h_residual, kBlack);

             // For simulation data, always use RMS
            // if (config.source == DataSource::SingleChipSim) {
            //     y_val = h_residual->GetRMS();
            //     y_err = h_residual->GetRMSError();
            // } else {
                if(fResidual) {
                    double chi2 = fResidual->GetChisquare();
                    int ndf = fResidual->GetNDF();
                    double chi2_per_ndf = (ndf > 0) ? (chi2 / ndf) : 9999.0;
                    y_rms = h_residual->GetRMS();

                    // if (chi2_per_ndf > 10) {
                    //     y_val = h_residual->GetRMS();
                    //     y_err = h_residual->GetRMSError();
                    // } else {
                    //     y_val = fResidual->GetParameter(2);
                    //     y_err = fResidual->GetParError(2);
                    // }

                    if(extract_tracking_resolution_) {
                        y_val = std::sqrt((y_rms*y_rms) - (tracking_resolution*tracking_resolution));
                    } else {
                        y_val = y_rms;
                    }
                    y_err = h_residual->GetRMSError();
                }
            //}
        } else if (quantity_to_extract == "clustersize") {
            y_val = h_clsize->GetMean();
            y_err = h_clsize->GetMeanError();
        }

        x_vals.push_back(x_val);
        y_vals.push_back(y_val);
        x_errs.push_back(0);
        y_errs.push_back(y_err);

        // inputFile->Close();
        // delete inputFile;
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
    //const std::pair<double, double>& x_range,
    const std::optional<std::pair<double, double>>& x_range) {
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

    if(configs[0].source == DataSource::SingleChipSim) {
        legend_line->Clear();
        legend_point->Clear();

        y_size = graphs.size() / 3 * 0.01;
        legend_y_min = 0.6;

        // TLegend* legend_point = new TLegend(0.15, legend_y_min, 0.90, 0.91);
        // legend_point->SetFillStyle(0);
        // legend_point->SetBorderSize(0);
        // legend_point->SetTextSize(0.03);
        legend_point->SetNColumns(3);

        // TLegend* legend_line = new TLegend(0.40, legend_y_min, 0.90, 0.91);
        // legend_line->SetFillStyle(0);
        // legend_line->SetBorderSize(0);
        // legend_line->SetTextSize(0.03);
        // legend_line->SetTextColor(kWhite);
        legend_line->SetNColumns(3);

        legend_point->SetX1NDC(0.15);
        legend_line->SetX1NDC(0.15);
        legend_point->SetY1NDC(legend_y_min);
        legend_line->SetY1NDC(legend_y_min);

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
    if(x_range) {
        mg->GetXaxis()->SetRangeUser(x_range->first, x_range->second);
    }
    //mg->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
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
        BEAM_INFO_ = "@KEK PF-AR Dec. 2024, 3 GeV/c electrons";
    } else if (config.source == DataSource::SPS202404) { // DataSource::SPS202404
        NAME_ = "sps202404";
        DUT_NAME_ = "CE65_6";
        BEAM_INFO_ = "@CERN SPS Apr. 2024, 120 GeV/c hadrons";
    } else if (config.source == DataSource::SingleChipSim) {
        NAME_ = "ce65sim202505";
        DUT_NAME_ = "CE65";
        BEAM_INFO_ = "3 GeV/c electrons";
    } else {
        LOG_ERROR.source("plot_BeamTest::check_residual_fits") << "No DataSource!";
        return;
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

            if(config.source == DataSource::SingleChipSim || config.source == DataSource::SingleChipDrift) {
                seed_thd_for_file = "0";
            }

            // //const std::string& thd = thresholds[current_index];
            // std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
            //                                    DATA_DIR_PATH_.c_str(),
            //                                    NAME_.c_str(),
            //                                    NAME_.c_str(),
            //                                    config.pixel_pitch.c_str(),
            //                                    config.chip_type.c_str(),
            //                                    config.voltage.c_str(),
            //                                    seed_thd_for_file.c_str(),  //config.seed_thd.c_str(),
            //                                    thd.c_str());
            
            // TFile* inputFile = TFile::Open(input_file_path.c_str());
            // if(!inputFile || inputFile->IsZombie()) {
            //     LOG_WARNING.source("plot_BeamTest::check_residual_fits") << "Cannot open file " << input_file_path;
            //     if(inputFile) delete inputFile;
            //     continue;
            // }

            // TH1D* h_res = (TH1D*)inputFile->Get(Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str()));
            // if(!h_res) {
            //     LOG_WARNING.source("plot_BeamTest::check_residual_fits") << "Histogram not found in " << input_file_path;
            //     delete inputFile;
            //     continue;
            // }

            // h_res->SetDirectory(0);
            // inputFile->Close();
            // delete inputFile;

            std::string base_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
                                                 DATA_DIR_PATH_.c_str(),
                                                 NAME_.c_str(),
                                                 NAME_.c_str(),
                                                 config.pixel_pitch.c_str(),
                                                 config.chip_type.c_str(),
                                                 config.voltage.c_str(),
                                                 seed_thd_for_file.c_str(),
                                                 thd.c_str());
            
            // if(config.source == DataSource::SingleChipSim) {

            // }
            
            std::string hist_name = Form("AnalysisCE65/%s/local_residuals/residualsX", DUT_NAME_.c_str());
            if(config.source == DataSource::SingleChipSim) {
                hist_name = "residual_x";
            }

            TH1D* h_res = get_merged_object<TH1D>(base_file_path, hist_name);
            
            if(!h_res) {
                LOG_WARNING.source("plot_BeamTest::check_residual_fits") << "Merged histogram not found for " << base_file_path;
                continue;
            }

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

std::vector<PlotConfig> plot_BeamTest::load_jsonConfigs(const nlohmann::json& j) {
    std::vector<PlotConfig> configs;

    try{
        for(const auto& item : j.at("comparison_configs")) {
            PlotConfig conf;
            std::string source_str = item.at("source");
            if(source_str == "KEK202412") conf.source = DataSource::KEK202412;
            else if(source_str == "SPS202404") conf.source = DataSource::SPS202404;
            else if(source_str == "SingleChipSim") conf.source = DataSource::SingleChipSim;
            else if(source_str == "SingleChipDrift") conf.source = DataSource::SingleChipDrift;
            else {
                LOG_WARNING.source("plot_BeamTest::load_jsonConfigs") << "Unknown data source in JSON: " << source_str;
                continue;
            }

            conf.pixel_pitch = item.at("pixel_pitch").get<std::string>();
            conf.chip_type = item.at("chip_type").get<std::string>();
            conf.voltage = item.at("voltage").get<std::string>();
            conf.seed_thd = item.at("seed_thd").get<std::string>();
            conf.scan_values = item.at("scan_values").get<std::vector<std::string>>();
            conf.legend_label = item.at("legend_label").get<std::string>();
            conf.marker_style = item.at("marker_style").get<int>();
            conf.marker_size = item.at("marker_size").get<double>();
            conf.line_style = item.at("line_style").get<int>();
            std::string color_string = item.at("color").get<std::string>();
            conf.color = string_to_ROOTColor(color_string);
            conf.adc_to_electron_factor = item.value("adc_to_electron_factor", 0.0);
            conf.base_file_name = item.value("base_file_name", "");
            conf.hist_path = item.value("hist_path", "");
            conf.hist_path_residual = item.value("hist_path_residual", "");
            conf.hist_path_clsize = item.value("hist_path_clsize", "");
            configs.push_back(conf);
        }
    } catch(nlohmann::json::exception& e) {
        LOG_ERROR.source("plot_BeamTest::load_jsonConfigs") << "JSON parse error: " << e.what();
    }
    return configs;
}

std::vector<InPixelPlotConfig> plot_BeamTest::load_jsonInPixelPlotConfigs(const nlohmann::json& j) {
    std::vector<InPixelPlotConfig> configs;
    // std::ifstream file(filename);
    // if(!file.is_open()) {
    //     LOG_ERROR.source("plot_BeamTest::load_jsonInPixelPlotConfigs") << "Could not open config file: " << filename;
    //     return configs;
    // }

    // nlohmann::json j;
    try{
        for(const auto& item : j.at("inpixel_plot_types")) {
            InPixelPlotConfig pt;
            pt.name = item.at("name").get<std::string>();
            pt.hist_path = item.at("hist_path").get<std::string>();
            pt.title = item.at("title").get<std::string>();
            pt.z_axis_title = item.at("z_axis_title").get<std::string>();
            pt.z_min = item.at("z_min").get<double>();
            pt.z_max = item.at("z_max").get<double>();
            if(item.contains("scale_factor")) {
                pt.scale_factor = item.at("scale_factor").get<double>();
            }
            configs.push_back(pt);
        }
    } catch (nlohmann::json::exception& e) {
            LOG_ERROR.source("plot_BeamTest::load_jsonInPixelPlotConfigs") << "JSON parse error: " << e.what();
    }
    return configs;
}

std::vector<PathConfig> plot_BeamTest::load_jsonPathConfigs(const nlohmann::json& j) {
    std::vector<PathConfig> configs;
    try {
        for (const auto& item : j.at("inpixel_path_configs")) {
            PathConfig pc;
            pc.name = item.at("name").get<std::string>();

            for(const auto& point_item : item.at("points")) {
                PathPoint p;
                p.label = point_item.at("label").get<std::string>();
                p.x_expr = point_item.at("x").get<std::string>();
                p.y_expr = point_item.at("y").get<std::string>();
                pc.points.push_back(p);
            }

            for(const auto& segment_item : item.at("segments")) { 
                PathSegment s;
                s.from = segment_item.at("from").get<std::string>();
                s.to = segment_item.at("to").get<std::string>();
                s.label = segment_item.at("label").get<std::string>();
                s.color_str = segment_item.at("color").get<std::string>();
                pc.segments.push_back(s);
            }
            configs.push_back(pc);
        }
    } catch (nlohmann::json::exception& e) {
        LOG_ERROR.source("plot_BeamTest::load_jsonPathConfigs") << "JSON parse error: " << e.what();
    }
    return configs;
}

void plot_BeamTest::drawBeamInfo(
    const std::string& beam_info,
    double x = 0.15,
    double y = 0.92
) {
    TLatex info;
    info.SetTextSize(0.035);
    info.SetTextFont(42);
    info.SetTextAlign(12);

    //info.DrawLatexNDC(x, y, "Beam Test")
    info.SetTextSize(0.03);
    info.DrawLatexNDC(x, y-0.04, beam_info.c_str());
    info.DrawLatexNDC(x, y-0.07, Form("Plotted on %s", TIME_.c_str()));
}

void plot_BeamTest::drawChipInfo(
    const ChipParameters& params,
    double x = 0.65,
    double y = 0.92
) {
    TLatex info;
    info.SetTextSize(0.028);
    info.SetTextFont(42);
    info.SetTextAlign(12);

    info.DrawLatexNDC(x, y, Form("Pitch: %s um", params.pixel_pitch.c_str()));
    info.DrawLatexNDC(x, y - 0.035, Form("Type: %s", params.chip_type.c_str()));
    info.DrawLatexNDC(x, y - 0.07, Form("Voltage: %s V", params.voltage.c_str()));
    info.DrawLatexNDC(x, y - 0.105, Form("Seed Thd: %s ADC", params.seed_thd.c_str()));
    info.DrawLatexNDC(x, y - 0.14, Form("Neighbor Thd: %s ADC", params.neighbor_thd.c_str()));
}

void plot_BeamTest::createAndSaveInPixelPlot(
    TH2D* hist,
    const InPixelPlotConfig& config,
    const ChipParameters& params,
    const std::string& beam_info,
    TDirectory* output_dir,
    const std::string& base_path,
    const std::string& chip_variation_name
) {
    canvas_->Clear();

    hist->Scale(config.scale_factor);
    hist->SetMinimum(config.z_min);
    hist->SetMaximum(config.z_max);
    hist->SetTitle(Form(";x w/in pixel [um];y w/in pixel [um];%s", config.z_axis_title.c_str()));

    canvas_->SetRightMargin(0.15);
    hist->Draw("COLZ");

    TLatex main_title;
    main_title.SetTextSize(0.045);
    main_title.SetTextFont(62);
    main_title.SetTextAlign(22);
    main_title.DrawLatexNDC(0.5, 0.96, config.title.c_str());

    drawBeamInfo(beam_info);
    drawChipInfo(params);

    std::string output_pdf_path = Form("%s/%s/%s_ce65_%s.pdf", 
                                       base_path.c_str(), 
                                       config.name.c_str(), 
                                       config.name.c_str(), 
                                       chip_variation_name.c_str());
    
    TString dir_path = gSystem->DirName(output_pdf_path.c_str());
    gSystem->mkdir(dir_path, kTRUE); // kTRUE makes it recursive

    output_dir->cd();
    //canvas_->Write(Form("%s_ce65_%s", config.name.c_str(), chip_variation_name.c_str()));
    canvas_->SaveAs(Form("%s/%s/%s_ce65_%s.pdf", base_path.c_str(), config.name.c_str(), config.name.c_str(), chip_variation_name.c_str()));
}

// void plot_BeamTest::run_inPixelAnalysis(
//     const std::vector<PlotConfig>& configs,
//     const std::vector<InPixelPlotConfig>& plot_types
// ) {
//     LOG_STATUS.source("plot_BeamTest::run_inPixelAnalysis") << "Start run for in-pixel analysis using JSON config.";

//     if(plot_types.empty()) {
//         LOG_ERROR.source("plot_BeamTest::run_inPixelAnalysis") << "No in-pixel plot types were difined in the JSON.";
//         return;
//     }

//     std::map<std::string, TFile*> output_files;
//     std::map<std::string, std::map<std::string, TDirectory*>> output_dirs;
//     gStyle->SetPalette(kViridis);

//     for(const auto& config : configs) {
//         if(config.source == DataSource::KEK202412) {
//             NAME_ = "kek202412";
//             DUT_NAME_ = "CE65_3";
//             BEAM_INFO_ = "e^{-} 3GeV/c @KEK PF-AR (Dec. 2024)";
//         } else if(config.source == DataSource::SPS202404) {
//             NAME_ = "sps202404";
//             DUT_NAME_ = "CE65_6";
//             BEAM_INFO_ = "hadron 120GeV/c @CERN SPS (Apr. 2024)";
//         } else {
//             LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Unknown data source enum value encountered.";
//             continue;
//         }

//         if(output_files.find(NAME_) == output_files.end()) {
//             std::string output_filename = "plot/experimentData_inPixel_" + NAME_ + ".root";
//             output_files[NAME_] = TFile::Open(output_filename.c_str(), "RECREATE");
//             for(const auto& pt : plot_types) {
//                 output_dirs[NAME_][pt.name] = output_files[NAME_]->mkdir(pt.name.c_str());
//             }
//         }

//         std::string data_dir_path = DATA_DIR_PATH_ + NAME_ + "/";
//         std::string plot_save_path = "plot/inPixel/" + NAME_;

//         for(const auto& neighbor_thd : config.scan_values) {
//             double seed_val = std::stod(config.seed_thd);
//             double neighbor_val = std::stod(neighbor_thd);
//             std::string seed_thd_for_file = (neighbor_val > seed_val) ? neighbor_thd : config.seed_thd;

//             ChipParameters current_params = {config.pixel_pitch, config.chip_type, config.voltage, seed_thd_for_file, neighbor_thd};

//             // std::string file_path = Form("%s%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
//             //                              data_dir_path.c_str(),
//             //                              NAME_.c_str(),
//             //                              config.pixel_pitch.c_str(),
//             //                              config.chip_type.c_str(),
//             //                              config.voltage.c_str(),
//             //                              seed_thd_for_file.c_str(),
//             //                              neighbor_thd.c_str());

//             // TFile* inputROOTFile = TFile::Open(file_path.c_str());
//             // if(!inputROOTFile || inputROOTFile->IsZombie()) {
//             //     LOG_ERROR.source("plot_BeamTest::run_inPixelAnalysis") << "Failed to open file: " << file_path;
//             //     if(inputROOTFile) delete inputROOTFile;
//             //     continue;
//             // }

//             std::string base_file_path = Form("%s%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
//                                                  data_dir_path.c_str(),
//                                                  NAME_.c_str(),
//                                                  config.pixel_pitch.c_str(),
//                                                  config.chip_type.c_str(),
//                                                  config.voltage.c_str(),
//                                                  seed_thd_for_file.c_str(),
//                                                  neighbor_thd.c_str());

//             std::string chip_variation_name = Form("%s_%s_%sV_SeedThd%se_NeighborThd%se",
//                                                    config.pixel_pitch.c_str(),
//                                                    config.chip_type.c_str(),
//                                                    config.voltage.c_str(),
//                                                    seed_thd_for_file.c_str(),
//                                                    neighbor_thd.c_str());
                                
//             for(auto plot_config : plot_types) {
//                 plot_config.hist_path = Form(plot_config.hist_path.c_str(), DUT_NAME_.c_str());

//                 // TProfile2D* prof = (TProfile2D*)inputROOTFile->Get(plot_config.hist_path.c_str());
//                 // if(!prof) {
//                 //     LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Histogram not found: " << plot_config.hist_path << "in file " << file_path;
//                 //     continue;
//                 // }

//                 // //TH2D* hist = prof->Projection2D();
//                 // TH2D* hist = plot_ExperimentData::convert_toTH2D(prof);
//                 // hist->SetDirectory(0);

//                 TProfile2D* prof = get_merged_object<TProfile2D>(base_file_path, plot_config.hist_path);
//                 if(!prof) {
//                     LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Merged TProfile2D not found: " << plot_config.hist_path << " for base " << base_file_path;
//                     continue;
//                 }

//                 TH2D* hist = plot_ExperimentData::convert_toTH2D(prof);
//                 hist->SetDirectory(0);
//                 delete prof;


//                 createAndSaveInPixelPlot(hist,
//                                          plot_config,
//                                          current_params,
//                                          BEAM_INFO_,
//                                          output_dirs[NAME_][plot_config.name],
//                                          plot_save_path,
//                                          chip_variation_name);

//                 delete hist;
//             }
//             // inputROOTFile->Close();
//             // delete inputROOTFile;
//         }
//     }
// }

void plot_BeamTest::run_inPixelAnalysis(
    const std::vector<PlotConfig>& configs,
    const std::vector<InPixelPlotConfig>& plot_types
) {
    LOG_STATUS.source("plot_BeamTest::run_inPixelAnalysis") << "Start run for in-pixel analysis using JSON config.";

    if(plot_types.empty()) {
        LOG_ERROR.source("plot_BeamTest::run_inPixelAnalysis") << "No in-pixel plot types were defined in the JSON.";
        return;
    }

    std::map<std::string, TFile*> output_files;
    std::map<std::string, std::map<std::string, TDirectory*>> output_dirs;
    gStyle->SetPalette(kViridis);

    for(const auto& config : configs) {
        if(config.source == DataSource::KEK202412) {
            NAME_ = "kek202412";
            DUT_NAME_ = "CE65_3";
            BEAM_INFO_ = "e^{-} 3GeV/c @KEK PF-AR (Dec. 2024)";
        } else if(config.source == DataSource::SPS202404) {
            NAME_ = "sps202404";
            DUT_NAME_ = "CE65_6";
            BEAM_INFO_ = "hadron 120GeV/c @CERN SPS (Apr. 2024)";
        } else {
            LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Unknown data source enum value encountered.";
            continue;
        }

        if(output_files.find(NAME_) == output_files.end()) {
            std::string output_filename = "plot/experimentData_inPixel_" + NAME_ + ".root";
            output_files[NAME_] = TFile::Open(output_filename.c_str(), "RECREATE");
            for(const auto& pt : plot_types) {
                output_dirs[NAME_][pt.name] = output_files[NAME_]->mkdir(pt.name.c_str());
            }
        }

        std::string data_dir_path = DATA_DIR_PATH_ + NAME_ + "/";
        std::string plot_save_path = "plot/inPixel/" + NAME_;

        for(const auto& neighbor_thd : config.scan_values) {
            double seed_val = std::stod(config.seed_thd);
            double neighbor_val = std::stod(neighbor_thd);
            std::string seed_thd_for_file = (neighbor_val > seed_val) ? neighbor_thd : config.seed_thd;

            ChipParameters current_params = {config.pixel_pitch, config.chip_type, config.voltage, seed_thd_for_file, neighbor_thd};

            // Use base_file_name from JSON if available
            std::string base_file_path;
            if (!config.base_file_name.empty()) {
                if (config.base_file_name.find("%s") != std::string::npos) {
                    base_file_path = Form(config.base_file_name.c_str(), seed_thd_for_file.c_str(), neighbor_thd.c_str());
                } else {
                    base_file_path = config.base_file_name;
                }
            } else {
                base_file_path = Form("%s%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
                                      data_dir_path.c_str(),
                                      NAME_.c_str(),
                                      config.pixel_pitch.c_str(),
                                      config.chip_type.c_str(),
                                      config.voltage.c_str(),
                                      seed_thd_for_file.c_str(),
                                      neighbor_thd.c_str());
            }

            std::string chip_variation_name = Form("%s_%s_%sV_SeedThd%se_NeighborThd%se",
                                                   config.pixel_pitch.c_str(),
                                                   config.chip_type.c_str(),
                                                   config.voltage.c_str(),
                                                   seed_thd_for_file.c_str(),
                                                   neighbor_thd.c_str());
                                        
            for(auto plot_config : plot_types) {
                plot_config.hist_path = Form(plot_config.hist_path.c_str(), DUT_NAME_.c_str());

                TProfile2D* prof = get_merged_object<TProfile2D>(base_file_path, plot_config.hist_path);
                if(!prof) {
                    LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Merged TProfile2D not found: " << plot_config.hist_path << " for base " << base_file_path;
                    continue;
                }

                TH2D* hist = plot_ExperimentData::convert_toTH2D(prof);
                hist->SetDirectory(0);
                delete prof;

                createAndSaveInPixelPlot(hist,
                                         plot_config,
                                         current_params,
                                         BEAM_INFO_,
                                         output_dirs[NAME_][plot_config.name],
                                         plot_save_path,
                                         chip_variation_name);

                delete hist;
            }
        }
    }
}

double plot_BeamTest::evaluate_expr(
    const std::string& expr,
    double half_pitch,
    double inset_x,
    double inset_y
) {
    if (expr == "0") return 0.0;
    if (expr == "half_pitch") return half_pitch;
    if (expr == "half_pitch/2") return half_pitch / 2.0;
    if (expr == "half_pitch - inset_x") return half_pitch - inset_x;
    if (expr == "half_pitch - inset_y") return half_pitch - inset_y;
    // 必要に応じて他の表現も追加
    try {
        return std::stod(expr);
    } catch (...) {
        LOG_ERROR << "Failed to evaluate expression: " << expr;
        return 0.0;
    }
    return 0;
}

// void plot_BeamTest::run_inPixelPathAnalysis(
//     const std::vector<PlotConfig>& plotConfigs,
//     const std::vector<InPixelPlotConfig>& inPixelPlotConfigs,
//     const std::vector<PathConfig>& pathConfigs
// ) {
//     LOG_STATUS.source("plot_BeamTest::run_inPixelPathAnalysis") << "Starting in-pixel path analysis.";
//     auto canvas_inpixel_ = std::make_unique<TCanvas>("canvas_inpixel_", "canvas_inpixel_", 2400, 800);

//     for (const auto& config : plotConfigs) {
//         if(config.source == DataSource::KEK202412) {
//             NAME_ = "kek202412";
//             DUT_NAME_ = "CE65_3";
//             BEAM_INFO_ = "e^{-} 3GeV/c @KEK PF-AR (Dec. 2024)";
//         } else if(config.source == DataSource::SPS202404) {
//             NAME_ = "sps202404";
//             DUT_NAME_ = "CE65_6";
//             BEAM_INFO_ = "hadron 120GeV/c @CERN SPS (Apr. 2024)";
//         }

//         LOG_STATUS.source("plot_BeamTest::run_inPixelPathAnalysis") << "Analysis is starting for " << NAME_;

//         for (const auto& neighbor_thd : config.scan_values) {
//             double seed_val = std::stod(config.seed_thd);
//             double neighbor_val = std::stod(neighbor_thd);
//             std::string seed_thd_for_file = (neighbor_val > seed_val) ? neighbor_thd : config.seed_thd;

//             std::string base_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
//                                               DATA_DIR_PATH_.c_str(),
//                                               NAME_.c_str(),
//                                               NAME_.c_str(),
//                                               config.pixel_pitch.c_str(),
//                                               config.chip_type.c_str(),
//                                               config.voltage.c_str(),
//                                               seed_thd_for_file.c_str(),
//                                               neighbor_thd.c_str());

//             for (const auto& plot_type : inPixelPlotConfigs) {
//                 std::string hist_path = Form(plot_type.hist_path.c_str(), DUT_NAME_.c_str());
//                 TProfile2D* prof = get_merged_object<TProfile2D>(base_file_path, hist_path);

//                 auto h_inpixel = std::unique_ptr<TH2D>(plot_ExperimentData::convert_toTH2D(prof));
//                 if (!h_inpixel) continue;

//                 for (const auto& path_config : pathConfigs) {
//                     LOG_INFO << "Processing path: " << path_config.name << " for plot: " << plot_type.name;

//                     canvas_inpixel_->Clear();
//                     TPad* pad_2d = new TPad("pad_map", "pad_map", 0.05, 0.0, 0.45, 1.0);
//                     pad_2d->SetMargin(0.10, 0.25, 0.12, 0.1);
//                     pad_2d->Draw();
//                     TPad* pad_1d = new TPad("pad_1d", "pad_1d", 0.45, 0.0, 0.95, 1.0);
//                     pad_1d->SetMargin(0.15, 0.2, 0.12, 0.1);
//                     pad_1d->SetGridy();
//                     pad_1d->Draw();

//                     pad_2d->cd();
//                     gStyle->SetPalette(kViridis);
//                     h_inpixel->SetStats(0);
//                     //h_inpixel->SetTitle(";In-pixel track intercept x [um];In-pixel track intercept y [um];#sqrt{#Delta x^{2} + #Delta y^{2}} (um)");
//                     h_inpixel->SetTitle(Form(";In-pixel track intercept x [um];In-pixel track intercept y [um];%s", plot_type.z_axis_title.c_str()));
//                     //h_inpixel->GetZaxis()->SetTitleOffset(1.2);
//                     //h_inpixel->GetYaxis()->SetLimits(3.5,9.5);
//                     h_inpixel->GetZaxis()->SetRangeUser(plot_type.z_min, plot_type.z_max);
//                     h_inpixel->GetXaxis()->SetTitleSize(0.06);
//                     h_inpixel->GetXaxis()->SetLabelSize(0.04);
//                     h_inpixel->GetXaxis()->SetTitleOffset(0.7);
//                     h_inpixel->GetYaxis()->SetTitleSize(0.06);
//                     h_inpixel->GetYaxis()->SetLabelSize(0.04);
//                     h_inpixel->GetYaxis()->SetTitleOffset(0.6);
//                     h_inpixel->GetXaxis()->SetTitleSize(0.06);
//                     h_inpixel->GetXaxis()->SetLabelSize(0.04);
//                     h_inpixel->GetXaxis()->SetTitleOffset(0.7);
//                     h_inpixel->GetZaxis()->SetTitleSize(0.06);
//                     h_inpixel->GetZaxis()->SetLabelSize(0.04);
//                     h_inpixel->GetZaxis()->SetTitleOffset(0.95);
//                     h_inpixel->GetXaxis()->SetTitleFont(42);
//                     h_inpixel->GetYaxis()->SetTitleFont(42);
//                     h_inpixel->GetZaxis()->SetTitleFont(42);
//                     h_inpixel->GetXaxis()->SetLabelFont(42);
//                     h_inpixel->GetYaxis()->SetLabelFont(42);
//                     h_inpixel->GetZaxis()->SetLabelFont(42);
//                     h_inpixel->Draw("COLZ");

//                     TPaletteAxis *palette = (TPaletteAxis*)h_inpixel->GetListOfFunctions()->FindObject("palette");
//                     if(palette) {
//                         double plot_right_edge = 1 - gPad->GetRightMargin();

//                         double gap = 0.025;
//                         double width = 0.04;

//                         palette->SetX1NDC(plot_right_edge + gap);
//                         palette->SetX2NDC(plot_right_edge + gap + width);

//                         gPad->Modified();
//                         gPad->Update();
//                     }
                    
//                     double pitch = std::stod(config.pixel_pitch);
//                     double half_pitch = pitch / 2.0;
//                     double inset_x = 0.1 * h_inpixel->GetXaxis()->GetBinWidth(1);
//                     double inset_y = 0.1 * h_inpixel->GetYaxis()->GetBinWidth(1);
                    
//                     std::map<std::string, TVector2> point_coords;
//                     for (const auto& p : path_config.points) {
//                         point_coords[p.label] = TVector2(
//                             evaluate_expr(p.x_expr, half_pitch, inset_x, inset_y),
//                             evaluate_expr(p.y_expr, half_pitch, inset_x, inset_y)
//                         );
//                     }
                    
//                     for (const auto& seg : path_config.segments) {
//                         TArrow arrow(point_coords.at(seg.from).X(), point_coords.at(seg.from).Y(),
//                                      point_coords.at(seg.to).X(), point_coords.at(seg.to).Y(), 0.015, ">");
//                         arrow.SetLineColorAlpha(string_to_ROOTColor(seg.color_str), 0.8);
//                         arrow.SetLineWidth(3);
//                         arrow.Draw();
//                     }
//                     for (const auto& p : path_config.points) {
//                         TMarker marker(point_coords.at(p.label).X(), point_coords.at(p.label).Y(), 24);
//                         marker.SetMarkerSize(3.5);
//                         // マーカーの色もJSONで定義可能にするとより良い
//                         marker.Draw("same");
//                         TLatex label_latex;
//                         label_latex.SetTextSize(0.06);
//                         label_latex.SetTextColor(kBlack);
//                         label_latex.SetTextAlign(22);
//                         label_latex.DrawLatex(point_coords.at(p.label).X(), point_coords.at(p.label).Y(), p.label);
//                     }

//                     // --- 1Dプロットのデータ抽出と描画 ---
//                     pad_1d->cd();
//                     TMultiGraph* mg_path = new TMultiGraph();
//                     mg_path->SetTitle(Form(";Distance along path [um];%s", plot_type.z_axis_title.c_str()));
//                     TLegend* legend_mg = new TLegend(0.8, 0.3, 1.0, 0.6);
//                     // ... 凡例のスタイリング ...
                    
//                     std::vector<std::unique_ptr<TGraphErrors>> graphs;
//                     double accumulated_distance = 0.0;
//                     std::vector<std::pair<double, std::string>> break_points;
                    
//                     // 開始点を記録
//                     if (!path_config.segments.empty()) {
//                          break_points.push_back({0.0, path_config.segments.front().from});
//                     }

//                     for (const auto& seg : path_config.segments) {
//                         const TVector2& p1 = point_coords.at(seg.from);
//                         const TVector2& p2 = point_coords.at(seg.to);
                        
//                         auto g_path = std::unique_ptr<TGraphErrors>(extract_data_along_path(h_inpixel.get(), p1.X(), p1.Y(), p2.X(), p2.Y(), accumulated_distance));
                        
//                         Color_t color = string_to_ROOTColor(seg.color_str);
//                         g_path->SetLineColor(color);
//                         g_path->SetLineWidth(2);

//                         auto g_err_band = std::unique_ptr<TGraphErrors>((TGraphErrors*)g_path->Clone());
//                         g_err_band->SetFillColorAlpha(color, 0.35);
//                         g_err_band->SetFillStyle(1001);
//                         mg_path->Add(g_err_band.get(), "3");
//                         mg_path->Add(g_path.get(), "L");

//                         legend_mg->AddEntry(g_path.get(), seg.label.c_str(), "lf");
//                         graphs.push_back(std::move(g_path));
                        
//                         accumulated_distance += (p2 - p1).Mod();
//                         break_points.push_back({accumulated_distance, seg.to});
//                     }

//                     mg_path->Draw("A");
//                     // ... 1Dプロットのスタイリング ...

//                     // 垂直線とラベルを描画
//                     TLine v_line; /* ... */ TLatex label; /* ... */
//                     // 重複を除いたユニークなブレークポイントを描画
//                     std::map<double, std::string> unique_breaks(break_points.begin(), break_points.end());
//                     for(const auto& bp : unique_breaks) {
//                         v_line.DrawLine(bp.first, plot_type.z_min, bp.first, plot_type.z_max);
//                         label.DrawLatex(bp.first, plot_type.z_max*1.05, bp.second.c_str());
//                     }
//                     legend_mg->Draw();
                    
//                     // ... ChipInfo描画 (変更なし) ...

//                     // 出力ファイル名にパス名を追加してユニークにする
//                     std::string output_filename = Form("plot/inPixel/%s/%s/%s_%s_inPixelPath_%s_%s_%s_%sV_Seed%s_Neighbor%s.pdf",
//                                                     NAME_.c_str(),
//                                                     plot_type.name.c_str(),
//                                                     path_config.name.c_str(), // ★追加
//                                                     NAME_.c_str(),
//                                                     plot_type.name.c_str(),
//                                                     config.pixel_pitch.c_str(),
//                                                     config.chip_type.c_str(),
//                                                     config.voltage.c_str(),
//                                                     seed_thd_for_file.c_str(),
//                                                     neighbor_thd.c_str());

//                     TString dir_path = gSystem->DirName(output_filename.c_str());
//                     gSystem->mkdir(dir_path, kTRUE);
//                     canvas_inpixel_->SaveAs(output_filename.c_str());
                    
//                     delete pad_1d;
//                     delete pad_2d;
//                     delete mg_path;
//                     delete legend_mg;
//                 } // end of path_configs loop
//             } // end of plot_types loop
//         } // end of scan_values loop
//     } // end of configs loop
// }

// ///////////////////////// not upgraded version
// void plot_BeamTest::run_inPixelPathAnalysis(const std::vector<PlotConfig>& configs, const std::vector<InPixelPlotConfig>& plot_types) {
//     LOG_STATUS.source("plot_BeamTest::run_inPixelPathAnalysis") << "Starting in-pixel path analysis.";
//     TCanvas* canvas_inpixel_ = new TCanvas("canvas_inpixel_", "canvas_inpixel_", 2400, 800);

//     for(const auto& config : configs) {
//         if(config.source == DataSource::KEK202412) {
//             NAME_ = "kek202412";
//             DUT_NAME_ = "CE65_3";
//             BEAM_INFO_ = "@KEK PF-AR Dec. 2024, 3 GeV/c electrons";
//         }else if(config.source == DataSource::SPS202404) {
//             NAME_ = "sps202404";
//             DUT_NAME_ = "CE65_6";
//             BEAM_INFO_ = "@CERN SPS Apr. 2024, 120 GeV/c hadrons";
//         }else if(config.source == DataSource::SingleChipSim) {
//             NAME_ = "ce65sim202505";
//             DUT_NAME_ = "CE65";
//             BEAM_INFO_ = "3 GeV/c electrons";
//         } else if(config.source == DataSource::SingleChipDrift) {
//             NAME_ = "ce65driftTime";
//             DUT_NAME_ = "CE65";
//             BEAM_INFO_ = "3 GeV/c electrons";
//         } else {
//             LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "No Datasource!";
//             return;
//         }

//         LOG_STATUS.source("plot_BeamTest::run_inPixelPathAnalysis") << "Analysis is starting for " << NAME_;

//         if(config.scan_values.empty()) {
//             LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "No scan_values found in config.";
//             return;
//         }
//         //const std::string& neighbor_thd = config.scan_values[0];
//         for(const auto& neighbor_thd : config.scan_values) {
//             double seed_val = std::stod(config.seed_thd);
//             double neighbor_val = std::stod(neighbor_thd);
//             std::string seed_thd_for_file = (neighbor_val > seed_val) ? neighbor_thd : config.seed_thd;

//             // std::string file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
//             //                             DATA_DIR_PATH_.c_str(),
//             //                             NAME_.c_str(),
//             //                             NAME_.c_str(),
//             //                             config.pixel_pitch.c_str(),
//             //                             config.chip_type.c_str(),
//             //                             config.voltage.c_str(),
//             //                             seed_thd_for_file.c_str(),
//             //                             neighbor_thd.c_str());
                                
//             // TFile* inputROOTFile = TFile::Open(file_path.c_str());
//             // if(!inputROOTFile || inputROOTFile->IsZombie()) {
//             //     LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "Failed to open file: " << file_path;
//             //     if(inputROOTFile) delete inputROOTFile;
//             //     return;
//             // }

//             if(config.source == DataSource::SingleChipSim || config.source == DataSource::SingleChipDrift) {
//                 seed_thd_for_file = "0";
//             }

//             std::string base_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
//                                               DATA_DIR_PATH_.c_str(),
//                                               NAME_.c_str(),
//                                               NAME_.c_str(),
//                                               config.pixel_pitch.c_str(),
//                                               config.chip_type.c_str(),
//                                               config.voltage.c_str(),
//                                               seed_thd_for_file.c_str(),
//                                               neighbor_thd.c_str());
            
//             // if(config.source == DataSource::SingleChipSim) {
//             //     base_file_path = Form("%s%s/n/ce65_p%s_%s_Thd%se_e3GeV_masetti",
//             //                           DATA_DIR_PATH_.c_str(),
//             //                           NAME_.c_str(),
//             //                           config.voltage.c_str(),
//             //                           config.pixel_pitch.c_str(),
//             //                           config.chip_type.c_str(),
//             //                           neighbor_thd.c_str());
//             // }

//             for(const auto plot_type : plot_types) {
//                 LOG_INFO.source("plot_BeamTest::run_inPixelPathAnalysis") << "Run for " << plot_type.name;

//                 std::string hist_path = Form(plot_type.hist_path.c_str(), DUT_NAME_.c_str());
                
//                 // TProfile2D* prof = (TProfile2D*)inputROOTFile->Get(hist_path.c_str());
//                 // if(!prof) {
//                 //     LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "Histogram not found: " << hist_path << " in " << file_path;
//                 //     inputROOTFile->Close();
//                 //     delete inputROOTFile;
//                 //     return;
//                 // }

//                 TProfile2D* prof = get_merged_object<TProfile2D>(base_file_path, hist_path);
//                 if(!prof) {
//                     LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Merged TProfile2D not found: " << hist_path << " for base " << base_file_path;
//                     continue;
//                 }

//                 TH2D* h_inpixel = plot_ExperimentData::convert_toTH2D(prof);
//                 h_inpixel->Scale(plot_type.scale_factor);
//                 h_inpixel->SetDirectory(0);
//                 delete prof;

//                 canvas_inpixel_->Clear();

//                 // TH2D* h_inpixel = plot_ExperimentData::convert_toTH2D(prof);
//                 // h_inpixel->Scale(plot_type.scale_factor);
//                 // h_inpixel->SetDirectory(0);
//                 // inputROOTFile->Close();
//                 // delete inputROOTFile;

//                 // TCanvas* canvas_inpixel_ = new TCanvas("canvas_inpixel_", "canvas_inpixel_", 1600, 700);

//                 TPad* pad_2d = new TPad("pad_map", "pad_map", 0.05, 0.0, 0.45, 1.0);
//                 pad_2d->SetMargin(0.10, 0.25, 0.12, 0.1);
//                 pad_2d->Draw();

//                 TPad* pad_1d = new TPad("pad_1d", "pad_1d", 0.45, 0.0, 0.95, 1.0);
//                 pad_1d->SetMargin(0.15, 0.2, 0.12, 0.1);
//                 pad_1d->SetGridy();
//                 pad_1d->Draw();

//                 pad_2d->cd();
//                 gStyle->SetPalette(kViridis);
//                 h_inpixel->SetStats(0);
//                 //h_inpixel->SetTitle(";In-pixel track intercept x [um];In-pixel track intercept y [um];#sqrt{#Delta x^{2} + #Delta y^{2}} (um)");
//                 h_inpixel->SetTitle(Form(";In-pixel track intercept x [um];In-pixel track intercept y [um];%s", plot_type.z_axis_title.c_str()));
//                 //h_inpixel->GetZaxis()->SetTitleOffset(1.2);
//                 //h_inpixel->GetYaxis()->SetLimits(3.5,9.5);
//                 h_inpixel->GetZaxis()->SetRangeUser(plot_type.z_min, plot_type.z_max);
//                 h_inpixel->GetXaxis()->SetTitleSize(0.06);
//                 h_inpixel->GetXaxis()->SetLabelSize(0.04);
//                 h_inpixel->GetXaxis()->SetTitleOffset(0.7);
//                 h_inpixel->GetYaxis()->SetTitleSize(0.06);
//                 h_inpixel->GetYaxis()->SetLabelSize(0.04);
//                 h_inpixel->GetYaxis()->SetTitleOffset(0.6);
//                 h_inpixel->GetXaxis()->SetTitleSize(0.06);
//                 h_inpixel->GetXaxis()->SetLabelSize(0.04);
//                 h_inpixel->GetXaxis()->SetTitleOffset(0.7);
//                 h_inpixel->GetZaxis()->SetTitleSize(0.06);
//                 h_inpixel->GetZaxis()->SetLabelSize(0.04);
//                 h_inpixel->GetZaxis()->SetTitleOffset(0.95);
//                 h_inpixel->GetXaxis()->SetTitleFont(42);
//                 h_inpixel->GetYaxis()->SetTitleFont(42);
//                 h_inpixel->GetZaxis()->SetTitleFont(42);
//                 h_inpixel->GetXaxis()->SetLabelFont(42);
//                 h_inpixel->GetYaxis()->SetLabelFont(42);
//                 h_inpixel->GetZaxis()->SetLabelFont(42);
//                 h_inpixel->Draw("COLZ");

//                 gPad->Update();
                
//                 TPaletteAxis *palette = (TPaletteAxis*)h_inpixel->GetListOfFunctions()->FindObject("palette");
//                 if(palette) {
//                     double plot_right_edge = 1 - gPad->GetRightMargin();

//                     double gap = 0.025;
//                     double width = 0.04;

//                     palette->SetX1NDC(plot_right_edge + gap);
//                     palette->SetX2NDC(plot_right_edge + gap + width);

//                     gPad->Modified();
//                     gPad->Update();
//                 }

//                 TLegend* box1 = new TLegend(0.12, 0.14, 0.60, 0.21);
//                 box1->SetFillColorAlpha(kWhite, 0.4);
//                 //box1->SetFillStyle(3001);
//                 box1->SetLineWidth(0);
//                 box1->Draw();

//                 TLatex info_latex;
//                 info_latex.SetTextFont(42);
//                 info_latex.SetTextSize(0.03);
//                 info_latex.DrawLatexNDC(0.13, 0.18, BEAM_INFO_.c_str());
//                 info_latex.DrawLatexNDC(0.13, 0.15, Form("Plotted on %s", TIME_.c_str()));

//                 double pitch = std::stod(config.pixel_pitch);
//                 double half_pitch = pitch / 2.0;

//                 //double inset = 0.00000001;
//                 double bin_width_x = h_inpixel->GetXaxis()->GetBinWidth(1);
//                 double bin_width_y = h_inpixel->GetYaxis()->GetBinWidth(1);

//                 double inset_x = 0.1 * bin_width_x;
//                 double inset_y = 0.1 * bin_width_y;

//                 double pA[2] = {0., 0.};
//                 double pB[2] = {0., half_pitch - inset_y};
//                 double pC[2] = {half_pitch - inset_x, half_pitch - inset_y};
//                 double pD[2] = {half_pitch - inset_x, 0.};

//                 // TLine* line_AB = new TLine(pA[0], pA[1], pB[0], pB[1]);
//                 // TLine* line_BC = new TLine(pB[0], pB[1], pC[0], pC[1]);
//                 // TLine* line_AD = new TLine(pA[0], pA[1], pD[0], pD[1]);
//                 // TLine* line_AC = new TLine(pA[0], pA[1], pC[0], pC[1]);
//                 double arrowSize = 0.015;
                
//                 TArrow* line_AB = new TArrow(pA[0], pA[1], pB[0], pB[1], arrowSize, ">");
//                 TArrow* line_BC = new TArrow(pB[0], pB[1], pC[0], pC[1], arrowSize, ">");
//                 TArrow* line_CA = new TArrow(pC[0], pC[1], pA[0], pA[1], arrowSize, ">");
//                 TArrow* line_AD = new TArrow(pA[0], pA[1], pD[0], pD[1], arrowSize, ">");

//                 TMarker* markerA = new TMarker(pA[0], pA[1], 24);
//                 TMarker* markerB = new TMarker(pB[0], pB[1], 24);
//                 TMarker* markerC = new TMarker(pC[0], pC[1], 24);
//                 TMarker* markerD = new TMarker(pD[0], pD[1], 24);

//                 std::vector<TMarker*> markers = {markerA, markerB, markerC, markerD};

//                 //std::vector<TLine*> lines = {line_AB, line_BC, line_AD, line_AC};
//                 std::vector<TArrow*> lines = {line_AB, line_BC, line_CA, line_AD};
//                 std::vector<Color_t> colors = {kOrange+7, kPink-3, kRed+1, kViolet-2};

//                 for(size_t i=0; i<markers.size(); ++i) {
//                     markers[i]->SetMarkerSize(3.5);
//                     markers[i]->SetMarkerColor(colors[i]);
//                     markers[i]->Draw("same");
//                 }

//                 for(size_t i=0; i<lines.size(); ++i) {
//                     lines[i]->SetLineColorAlpha(colors[i], 0.8);
//                     lines[i]->SetLineWidth(3);
//                     lines[i]->SetLineStyle(1);
//                     lines[i]->Draw();
//                 }

//                 TLatex label_latex;
//                 label_latex.SetTextSize(0.06);
//                 label_latex.SetTextColor(kBlack);
//                 label_latex.SetTextAlign(22);
//                 label_latex.DrawLatex(pA[0], pA[1] - 1.0, "A");
//                 label_latex.DrawLatex(pB[0], pB[1] + 1.0, "B");
//                 label_latex.DrawLatex(pC[0], pC[1] + 1.0, "C");
//                 label_latex.DrawLatex(pD[0], pD[1] - 1.0, "D");

//                 auto extract_data_along_path = 
//                     [&](TH2D* hist, double x1, double y1, double x2, double y2, double dist_offset) ->TGraphErrors* {
//                     std::vector<double> dist, val, err_dist, err_val;

//                     TAxis* xAxis = hist->GetXaxis();
//                     TAxis* yAxis = hist->GetYaxis();

//                     int n_steps = 100;
//                     for(int i=0; i<=n_steps; ++i) {
//                         double t = (double)i / n_steps;
//                         double current_x = x1 + t * (x2 - x1);
//                         double current_y = y1 + t * (y2 - y1);

//                         int bin = hist->FindBin(current_x, current_y);

//                         if(i > 0 && hist->FindBin(x1 + (double)(i-1)/n_steps * (x2-x1), y1 + (double)(i-1)/n_steps * (y2-y1)) == bin) {
//                             continue;
//                         }

//                         double current_dist = dist_offset + std::sqrt(pow(current_x - x1, 2) + pow(current_y - y1, 2));

//                         dist.push_back(current_dist);
//                         val.push_back(hist->GetBinContent(bin));
//                         err_dist.push_back(0);
//                         err_val.push_back(hist->GetBinError(bin));
//                     }

//                     if(dist.empty()) return new TGraphErrors();
//                     return new TGraphErrors(dist.size(), dist.data(), val.data(), err_dist.data(), err_val.data());
//                 };

//                 double dist_AB = std::sqrt(pow(pB[0] - pA[0], 2) + pow(pB[1] - pA[1], 2));
//                 double dist_BC = std::sqrt(pow(pB[0] - pC[0], 2) + pow(pB[1] - pC[1], 2));
//                 double dist_AC = std::sqrt(pow(pA[0] - pC[0], 2) + pow(pA[1] - pC[1], 2));
//                 double dist_AD = std::sqrt(pow(pA[0] - pD[0], 2) + pow(pA[1] - pD[1], 2));
//                 // A-B-C-A-D
//                 TGraphErrors* g_path_AB = extract_data_along_path(h_inpixel, pA[0], pA[1], pB[0], pB[1], 0.0);
//                 TGraphErrors* g_path_BC = extract_data_along_path(h_inpixel, pB[0], pB[1], pC[0], pC[1], dist_AB);
//                 TGraphErrors* g_path_CA = extract_data_along_path(h_inpixel, pC[0], pC[1], pA[0], pA[1], dist_AB + dist_BC);
//                 TGraphErrors* g_path_AD = extract_data_along_path(h_inpixel, pA[0], pA[1], pD[0], pD[1], dist_AB + dist_AC + dist_BC);

//                 pad_1d->cd();
//                 TMultiGraph* mg_path = new TMultiGraph();
//                 //mg_path->SetTitle(Form(";Distance along the path [um];#sqrt{#Delta x^{2} + #Delta y^{2}} [um]"));
//                 mg_path->SetTitle(Form(";Distance along the path [um];%s", plot_type.z_axis_title.c_str()));
                
//                 std::vector<double> x_coords = {
//                     0.0, // A
//                     dist_AB, // B
//                     dist_AB + dist_BC, // C
//                     dist_AB + dist_BC + dist_AC, // A (return)
//                     dist_AB + dist_BC + dist_AC + dist_AD // D
//                 };

//                 std::vector<std::string> point_labels = {"A", "B", "C", "A", "D"};
//                 TLine v_line;
//                 v_line.SetLineStyle(2);
//                 v_line.SetLineColor(kGray+2);
//                 v_line.SetLineWidth(1);

//                 TLatex label;
//                 label.SetTextSize(0.04);
//                 label.SetTextAlign(23);

//                 // for(size_t i=0; i<x_coords.size(); i++) {
//                 //     v_line.DrawLine(x_coords[i], y_min_range, x_coords[i], y_max_range);
//                 //     label.DrawLatex(x_coords[i], y_max_range, point_labels[i].c_str());
//                 // }

//                 std::vector<TGraphErrors*> graphs = {g_path_AB, g_path_BC, g_path_CA, g_path_AD};
//                 std::vector<std::string> labels = {"A #rightarrow B", "B #rightarrow C", "C #rightarrow A", "A #rightarrow D"};

//                 for(size_t i=0; i<graphs.size(); ++i) {
//                     graphs[i]->SetLineColor(colors[i]);
//                     graphs[i]->SetLineWidth(2);

//                     TGraphErrors* g_err_band = (TGraphErrors*)graphs[i]->Clone();
//                     g_err_band->SetFillColorAlpha(colors[i], 0.35);
//                     g_err_band->SetFillStyle(1001); // Solid fill
//                     mg_path->Add(g_err_band, "3");

//                     mg_path->Add(graphs[i], "L");
//                 }

//                 mg_path->Draw("A");
//                 // mg_path->GetYaxis()->SetRangeUser(0, 9);
//                 // mg_path->GetXaxis()->SetTitleOffset(1.1);
//                 // mg_path->GetYaxis()->SetTitleOffset(1.3);
//                 //mg_path->GetYaxis()->SetLimits(y_min_range, y_max_range);
//                 mg_path->GetYaxis()->SetRangeUser(plot_type.z_min, plot_type.z_max);
//                 mg_path->GetXaxis()->SetTitleSize(0.06);
//                 mg_path->GetXaxis()->SetLabelSize(0.04);
//                 mg_path->GetXaxis()->SetTitleOffset(0.7);
//                 mg_path->GetYaxis()->SetTitleSize(0.06);
//                 mg_path->GetYaxis()->SetLabelSize(0.04);
//                 mg_path->GetYaxis()->SetTitleOffset(0.8);
//                 mg_path->GetXaxis()->SetTitleFont(42);
//                 mg_path->GetYaxis()->SetTitleFont(42);
//                 mg_path->GetXaxis()->SetLabelFont(42);
//                 mg_path->GetYaxis()->SetLabelFont(42);

//                 TLegend* legend_mg = new TLegend(0.8, 0.3, 1.0, 0.6);
//                 legend_mg->SetFillStyle(0);
//                 legend_mg->SetTextSize(0.04);
//                 legend_mg->SetBorderSize(0);
//                 for(size_t i=0; i<labels.size(); ++i) {
//                     legend_mg->AddEntry(graphs[i], labels[i].c_str(), "lf");
//                 }
//                 legend_mg->Draw();

//                 // double y_line_min = mg_path->GetYaxis()->GetXmin();
//                 // double y_line_max = mg_path->GetYaxis()->GetXmax();
//                 for(size_t i=0; i<x_coords.size(); i++) {
//                     v_line.DrawLine(x_coords[i], plot_type.z_min, x_coords[i], plot_type.z_max);
//                     label.DrawLatex(x_coords[i], plot_type.z_max*1.05, point_labels[i].c_str());
//                 }

//                 ChipParameters params = {config.pixel_pitch, config.chip_type, config.voltage, seed_thd_for_file, neighbor_thd};
//                 drawChipInfo(params, 0.81, 0.88);

//                 std::string safe_legend = config.legend_label;
//                 std::replace(safe_legend.begin(), safe_legend.end(), ' ', '_');
//                 std::replace(safe_legend.begin(), safe_legend.end(), ',', '_'); 

//                 std::string output_filename = Form("plot/inPixel/%s/%s/%s_%s_inPixelPath_%s_%s_%sV_Seed%s_Neighbor%s.pdf",
//                                                 NAME_.c_str(),
//                                                 plot_type.name.c_str(),
//                                                 NAME_.c_str(),
//                                                 plot_type.name.c_str(),
//                                                 config.pixel_pitch.c_str(),
//                                                 config.chip_type.c_str(),
//                                                 config.voltage.c_str(),
//                                                 seed_thd_for_file.c_str(),
//                                                 neighbor_thd.c_str());

//                 TString dir_path = gSystem->DirName(output_filename.c_str());
//                 gSystem->mkdir(dir_path, kTRUE);
//                 canvas_inpixel_->SaveAs(output_filename.c_str());
//             }

//             // inputROOTFile->Close();
//             // delete inputROOTFile;
//         }
//     }
// }

void plot_BeamTest::run_inPixelPathAnalysis(const std::vector<PlotConfig>& configs, const std::vector<InPixelPlotConfig>& plot_types) {
    LOG_STATUS.source("plot_BeamTest::run_inPixelPathAnalysis") << "Starting in-pixel path analysis.";
    TCanvas* canvas_inpixel_ = new TCanvas("canvas_inpixel_", "canvas_inpixel_", 2400, 800);

    for(const auto& config : configs) {
        if(config.source == DataSource::KEK202412) {
            NAME_ = "kek202412";
            DUT_NAME_ = "CE65_3";
            BEAM_INFO_ = "@KEK PF-AR Dec. 2024, 3 GeV/c electrons";
        }else if(config.source == DataSource::SPS202404) {
            NAME_ = "sps202404";
            DUT_NAME_ = "CE65_6";
            BEAM_INFO_ = "@CERN SPS Apr. 2024, 120 GeV/c hadrons";
        }else if(config.source == DataSource::SingleChipSim) {
            NAME_ = "ce65inPixel2601";
            DUT_NAME_ = "CE65";
            BEAM_INFO_ = "120 GeV/c pions";
        } else if(config.source == DataSource::SingleChipDrift) {
            NAME_ = "ce65driftTime";
            DUT_NAME_ = "CE65";
            BEAM_INFO_ = "120 GeV/c pions";
        } else {
            LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "No Datasource!";
            return;
        }

        LOG_STATUS.source("plot_BeamTest::run_inPixelPathAnalysis") << "Analysis is starting for " << NAME_;

        if(config.scan_values.empty()) {
            LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "No scan_values found in config.";
            return;
        }

        for(const auto& neighbor_thd : config.scan_values) {
            double seed_val = std::stod(config.seed_thd);
            double neighbor_val = std::stod(neighbor_thd);
            std::string seed_thd_for_file = (neighbor_val > seed_val) ? neighbor_thd : config.seed_thd;

            if(config.source == DataSource::SingleChipSim || config.source == DataSource::SingleChipDrift) {
                //seed_thd_for_file = "0";
            }

            // Use base_file_name from JSON if available
            std::string base_file_path;
            if (!config.base_file_name.empty()) {
                if (config.base_file_name.find("%s") != std::string::npos) {
                    base_file_path = Form(config.base_file_name.c_str(), seed_thd_for_file.c_str(), neighbor_thd.c_str());
                } else {
                    base_file_path = config.base_file_name;
                }
            } else {
                base_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se",
                                      DATA_DIR_PATH_.c_str(),
                                      NAME_.c_str(),
                                      NAME_.c_str(),
                                      config.pixel_pitch.c_str(),
                                      config.chip_type.c_str(),
                                      config.voltage.c_str(),
                                      seed_thd_for_file.c_str(),
                                      neighbor_thd.c_str());
            }
            
            for(const auto plot_type : plot_types) {
                LOG_INFO.source("plot_BeamTest::run_inPixelPathAnalysis") << "Run for " << plot_type.name;

                std::string hist_path = Form(plot_type.hist_path.c_str(), DUT_NAME_.c_str());
                
                TProfile2D* prof = get_merged_object<TProfile2D>(base_file_path, hist_path);
                if(!prof) {
                    LOG_WARNING.source("plot_BeamTest::run_inPixelAnalysis") << "Merged TProfile2D not found: " << hist_path << " for base " << base_file_path;
                    continue;
                }

                TH2D* h_inpixel = plot_ExperimentData::convert_toTH2D(prof);
                h_inpixel->Scale(plot_type.scale_factor);
                h_inpixel->SetDirectory(0);
                delete prof;

                h_inpixel->Rebin2D(2, 2);
                h_inpixel->Scale(1.0 / 4);

                canvas_inpixel_->Clear();

                TPad* pad_2d = new TPad("pad_map", "pad_map", 0.05, 0.0, 0.45, 1.0);
                pad_2d->SetMargin(0.10, 0.25, 0.12, 0.1);
                pad_2d->Draw();

                TPad* pad_1d = new TPad("pad_1d", "pad_1d", 0.45, 0.0, 0.95, 1.0);
                pad_1d->SetMargin(0.15, 0.2, 0.12, 0.1);
                pad_1d->SetGridy();
                pad_1d->Draw();

                pad_2d->cd();
                gStyle->SetPalette(kViridis);
                h_inpixel->SetStats(0);
                h_inpixel->SetTitle(Form(";In-pixel track intercept x [um];In-pixel track intercept y [um];%s", plot_type.z_axis_title.c_str()));
                h_inpixel->GetZaxis()->SetRangeUser(plot_type.z_min, plot_type.z_max);
                h_inpixel->GetXaxis()->SetTitleSize(0.06);
                h_inpixel->GetXaxis()->SetLabelSize(0.04);
                h_inpixel->GetXaxis()->SetTitleOffset(0.7);
                h_inpixel->GetYaxis()->SetTitleSize(0.06);
                h_inpixel->GetYaxis()->SetLabelSize(0.04);
                h_inpixel->GetYaxis()->SetTitleOffset(0.6);
                h_inpixel->GetZaxis()->SetTitleSize(0.06);
                h_inpixel->GetZaxis()->SetLabelSize(0.04);
                h_inpixel->GetZaxis()->SetTitleOffset(0.95);
                h_inpixel->Draw("COLZ");

                gPad->Update();
                
                TPaletteAxis *palette = (TPaletteAxis*)h_inpixel->GetListOfFunctions()->FindObject("palette");
                if(palette) {
                    double plot_right_edge = 1 - gPad->GetRightMargin();
                    double gap = 0.025;
                    double width = 0.04;
                    palette->SetX1NDC(plot_right_edge + gap);
                    palette->SetX2NDC(plot_right_edge + gap + width);
                    gPad->Modified();
                    gPad->Update();
                }

                TLegend* box1 = new TLegend(0.12, 0.14, 0.60, 0.21);
                box1->SetFillColorAlpha(kWhite, 0.4);
                box1->SetLineWidth(0);
                box1->Draw();

                TLatex info_latex;
                info_latex.SetTextFont(42);
                info_latex.SetTextSize(0.03);
                info_latex.DrawLatexNDC(0.13, 0.18, BEAM_INFO_.c_str());
                info_latex.DrawLatexNDC(0.13, 0.15, Form("Plotted on %s", TIME_.c_str()));

                double pitch = std::stod(config.pixel_pitch);
                double half_pitch = pitch / 2.0;
                double bin_width_x = h_inpixel->GetXaxis()->GetBinWidth(1);
                double bin_width_y = h_inpixel->GetYaxis()->GetBinWidth(1);
                double inset_x = 0.1 * bin_width_x;
                double inset_y = 0.1 * bin_width_y;

                double pA[2] = {0., 0.};
                double pB[2] = {0., half_pitch - inset_y};
                double pC[2] = {half_pitch - inset_x, half_pitch - inset_y};
                double pD[2] = {half_pitch - inset_x, 0.};

                double arrowSize = 0.015;
                TArrow* line_AB = new TArrow(pA[0], pA[1], pB[0], pB[1], arrowSize, ">");
                TArrow* line_BC = new TArrow(pB[0], pB[1], pC[0], pC[1], arrowSize, ">");
                TArrow* line_CA = new TArrow(pC[0], pC[1], pA[0], pA[1], arrowSize, ">");
                TArrow* line_AD = new TArrow(pA[0], pA[1], pD[0], pD[1], arrowSize, ">");

                TMarker* markerA = new TMarker(pA[0], pA[1], 24);
                TMarker* markerB = new TMarker(pB[0], pB[1], 24);
                TMarker* markerC = new TMarker(pC[0], pC[1], 24);
                TMarker* markerD = new TMarker(pD[0], pD[1], 24);

                std::vector<TMarker*> markers = {markerA, markerB, markerC, markerD};
                std::vector<TArrow*> lines = {line_AB, line_BC, line_CA, line_AD};
                std::vector<Color_t> colors = {kOrange+7, kPink-3, kRed+1, kViolet-2};

                for(size_t i=0; i<markers.size(); ++i) {
                    markers[i]->SetMarkerSize(3.5);
                    markers[i]->SetMarkerColor(colors[i]);
                    markers[i]->Draw("same");
                }

                for(size_t i=0; i<lines.size(); ++i) {
                    lines[i]->SetLineColorAlpha(colors[i], 0.8);
                    lines[i]->SetLineWidth(3);
                    lines[i]->Draw();
                }

                TLatex label_latex;
                label_latex.SetTextSize(0.06);
                label_latex.SetTextAlign(22);
                label_latex.DrawLatex(pA[0], pA[1] - 1.0, "A");
                label_latex.DrawLatex(pB[0], pB[1] + 1.0, "B");
                label_latex.DrawLatex(pC[0], pC[1] + 1.0, "C");
                label_latex.DrawLatex(pD[0], pD[1] - 1.0, "D");

                auto extract_data_along_path = 
                    [&](TH2D* hist, double x1, double y1, double x2, double y2, double dist_offset) ->TGraphErrors* {
                    std::vector<double> dist, val, err_dist, err_val;
                    int n_steps = 100;
                    for(int i=0; i<=n_steps; ++i) {
                        double t = (double)i / n_steps;
                        double current_x = x1 + t * (x2 - x1);
                        double current_y = y1 + t * (y2 - y1);
                        int bin = hist->FindBin(current_x, current_y);
                        if(i > 0 && hist->FindBin(x1 + (double)(i-1)/n_steps * (x2-x1), y1 + (double)(i-1)/n_steps * (y2-y1)) == bin) continue;
                        double current_dist = dist_offset + std::sqrt(pow(current_x - x1, 2) + pow(current_y - y1, 2));
                        dist.push_back(current_dist);
                        val.push_back(hist->GetBinContent(bin));
                        err_dist.push_back(0);
                        err_val.push_back(hist->GetBinError(bin));
                    }
                    return new TGraphErrors(dist.size(), dist.data(), val.data(), err_dist.data(), err_val.data());
                };

                double dist_AB = std::sqrt(pow(pB[0] - pA[0], 2) + pow(pB[1] - pA[1], 2));
                double dist_BC = std::sqrt(pow(pB[0] - pC[0], 2) + pow(pB[1] - pC[1], 2));
                double dist_AC = std::sqrt(pow(pA[0] - pC[0], 2) + pow(pA[1] - pC[1], 2));
                double dist_AD = std::sqrt(pow(pA[0] - pD[0], 2) + pow(pA[1] - pD[1], 2));

                TGraphErrors* g_path_AB = extract_data_along_path(h_inpixel, pA[0], pA[1], pB[0], pB[1], 0.0);
                TGraphErrors* g_path_BC = extract_data_along_path(h_inpixel, pB[0], pB[1], pC[0], pC[1], dist_AB);
                TGraphErrors* g_path_CA = extract_data_along_path(h_inpixel, pC[0], pC[1], pA[0], pA[1], dist_AB + dist_BC);
                TGraphErrors* g_path_AD = extract_data_along_path(h_inpixel, pA[0], pA[1], pD[0], pD[1], dist_AB + dist_AC + dist_BC);

                pad_1d->cd();
                TMultiGraph* mg_path = new TMultiGraph();
                mg_path->SetTitle(Form(";Distance along the path [um];%s", plot_type.z_axis_title.c_str()));
                
                std::vector<double> x_coords = {0.0, dist_AB, dist_AB + dist_BC, dist_AB + dist_BC + dist_AC, dist_AB + dist_BC + dist_AC + dist_AD};
                std::vector<std::string> point_labels = {"A", "B", "C", "A", "D"};
                
                std::vector<TGraphErrors*> graphs = {g_path_AB, g_path_BC, g_path_CA, g_path_AD};
                std::vector<std::string> labels = {"A #rightarrow B", "B #rightarrow C", "C #rightarrow A", "A #rightarrow D"};

                for(size_t i=0; i<graphs.size(); ++i) {
                    graphs[i]->SetLineColor(colors[i]);
                    graphs[i]->SetLineWidth(2);
                    TGraphErrors* g_err_band = (TGraphErrors*)graphs[i]->Clone();
                    g_err_band->SetFillColorAlpha(colors[i], 0.35);
                    g_err_band->SetFillStyle(1001);
                    mg_path->Add(g_err_band, "3");
                    mg_path->Add(graphs[i], "L");
                }

                mg_path->Draw("A");
                mg_path->GetYaxis()->SetRangeUser(plot_type.z_min, plot_type.z_max);
                mg_path->GetXaxis()->SetTitleSize(0.06);
                mg_path->GetYaxis()->SetTitleSize(0.06);

                TLegend* legend_mg = new TLegend(0.8, 0.3, 1.0, 0.6);
                legend_mg->SetFillStyle(0);
                legend_mg->SetBorderSize(0);
                for(size_t i=0; i<labels.size(); ++i) legend_mg->AddEntry(graphs[i], labels[i].c_str(), "lf");
                legend_mg->Draw();

                TLine v_line;
                v_line.SetLineStyle(2);
                v_line.SetLineColor(kGray+2);
                TLatex label;
                label.SetTextSize(0.04);
                label.SetTextAlign(23);
                for(size_t i=0; i<x_coords.size(); i++) {
                    v_line.DrawLine(x_coords[i], plot_type.z_min, x_coords[i], plot_type.z_max);
                    label.DrawLatex(x_coords[i], plot_type.z_max*1.05, point_labels[i].c_str());
                }

                ChipParameters params = {config.pixel_pitch, config.chip_type, config.voltage, seed_thd_for_file, neighbor_thd};
                drawChipInfo(params, 0.81, 0.88);

                std::string output_filename = Form("plot/inPixel/%s/%s/%s_%s_inPixelPath_%s_%s_%sV_Seed%s_Neighbor%s.pdf",
                                                NAME_.c_str(), plot_type.name.c_str(), NAME_.c_str(), plot_type.name.c_str(),
                                                config.pixel_pitch.c_str(), config.chip_type.c_str(), config.voltage.c_str(),
                                                seed_thd_for_file.c_str(), neighbor_thd.c_str());

                TString dir_path = gSystem->DirName(output_filename.c_str());
                gSystem->mkdir(dir_path, kTRUE);
                canvas_inpixel_->SaveAs(output_filename.c_str());

                // Output file path
                std::string root_output_name = Form("plot/inPixelPathAnalysis_%s.root", NAME_.c_str());
                TFile* outFile = TFile::Open(root_output_name.c_str(), "UPDATE");

                if (outFile && !outFile->IsZombie()) {
                    // 1. Create or get Plot Type directory (e.g., cluster_size, residual)
                    TDirectory* typeDir = outFile->mkdir(plot_type.name.c_str(), "", true);
                    typeDir->cd();

                    // 2. Create or get Threshold condition directory
                    std::string sub_dir_name = Form("Seed%s_Neighbor%s", seed_thd_for_file.c_str(), neighbor_thd.c_str());
                    TDirectory* condDir = typeDir->mkdir(sub_dir_name.c_str(), "", true);
                    condDir->cd();

                    // 3. Define a short name for the quantity to include in object names
                    std::string q_name = plot_type.name;

                    // 4. Save the 2D Map
                    // Use kOverwrite to replace the existing object instead of creating a new cycle (;2, ;3, etc.)
                    h_inpixel->Write(Form("map_2d_%s", q_name.c_str()), TObject::kOverwrite);

                    // 5. Save the 1D Path Graphs with explicit path identifiers
                    g_path_AB->Write(Form("path_A_to_B_%s", q_name.c_str()), TObject::kOverwrite);
                    g_path_BC->Write(Form("path_B_to_C_%s", q_name.c_str()), TObject::kOverwrite);
                    g_path_CA->Write(Form("path_C_to_A_%s", q_name.c_str()), TObject::kOverwrite);
                    g_path_AD->Write(Form("path_A_to_D_%s", q_name.c_str()), TObject::kOverwrite);

                    // 6. Save the MultiGraph
                    mg_path->Write(Form("mg_combined_path_%s", q_name.c_str()), TObject::kOverwrite);

                    outFile->Close();
                    delete outFile;
                } else {
                    LOG_ERROR.source("plot_BeamTest::run_inPixelPathAnalysis") << "Could not open ROOT file for writing: " << root_output_name;
                }
            }
        }
    }
}

void plot_BeamTest::BeamTest_main(int argc, char* argv[]) {
    cxxopts::Options options("plot_BeamTest::BeamTest_main", "Beamtest object code");
    options.add_options()
        ("f,file", "Json file name", cxxopts::value<std::string>())
        ("i,inpixel", "Run in-Pixel analysis instead of comparison plots.")
        ("e,electrons", "Scale threshold axis to electrons (e-)")
        ("t,tracking_resolution", "residual tracking resolution for truly position resolution")
        ("h,help", "show help message");
    
    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        std::cout << options.help() << std::endl;
        return;
    }

    if (result.count("electrons")) {
        use_electron_scale_ = true;
        LOG_STATUS.source("plot_BeamTest::BeamTest_main") << "Threshold scaling to [e-] enabled.";
    }
    if(result.count("tracking_resolution")) {
        extract_tracking_resolution_ = true;
        LOG_STATUS.source("plot_BeamTest::BeamTest_main") << "Extract tracking resolution from position resolution.";
    }

    std::string config_filename;
    if(result.count("file")) {
        config_filename = result["file"].as<std::string>();
    } else {
        LOG_WARNING.source("plot_BeamTest::BeamTest_main") << "JSON file name is NOT defined. Using default file name: plot_BeamTest/json/analysis.json";
        config_filename = "plot_BeamTest/json/analysis.json";
    }

    std::ifstream file(config_filename);
    nlohmann::json j;
    file >> j;

    if(j.is_null() || !j.is_object()) {
        LOG_ERROR.source("plot_BeamTest::BeamTest_main") << "Failed to parse JSON or file is empty or invalid: " << config_filename;
        return;
    }

    // std::vector<PlotConfig> configs = plot_BeamTest::load_jsonConfigs(config_filename);

    // if(configs.empty()) {
    //     LOG_ERROR.source("plot_BeamTest::BeamTest_main") << "No configurations were loaded.";
    //     return;
    // }

    if(result.count("inpixel")) {
        LOG_STATUS.source("plot_BeamTest::BeamTest_main") << "Executing In-Pixel Analysis based on " << config_filename;
        std::vector<PlotConfig> configs = load_jsonConfigs(j);
        std::vector<InPixelPlotConfig> inpixel_configs = load_jsonInPixelPlotConfigs(j);
        
        //run_inPixelAnalysis(configs, inpixel_configs);
        run_inPixelPathAnalysis(configs, inpixel_configs);
        LOG_STATUS.source("plot_BeamTest::BeamTest_main") << "In-Pixel Analysis finished.";
    } else {
        LOG_STATUS.source("plot_BeamTest::BeamTest_main") << "Executing Comparison Plots based on " << config_filename;
        std::vector<PlotConfig> configs = load_jsonConfigs(j);
        
        run_plots(configs);
        LOG_STATUS.source("Comparison Plots finished.");
    }

    //plot_BeamTest::run_plots(my_comparison);
}