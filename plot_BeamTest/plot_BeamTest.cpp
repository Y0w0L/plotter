#include "plot_BeamTest.h"

// Constructer
plot_BeamTest::plot_BeamTest(DataSource source) : source_(source) {
    LOG_STATUS.source("plot_BeamTest::plot_BeamTest") << "plot_BeamTest object is created.";

    if (source == DataSource::KEK202412) {
        NAME_ = "kek202412";
        DUT_NAME_ = "CE65_3";
        BEAM_INFO_ = "e^{-} 3GeV/c @KEK-PFAR (Dec. 2024)";
    } else if (source == DataSource::SPS202404) {
        NAME_ = "sps202404";
        DUT_NAME_ = "CE65_6";
        BEAM_INFO_ = "hadron 120GeV/c @CERN-SPS (Apr. 2024)";
    }
    else {
        LOG_ERROR.source("plot_BeamTest::plot_BeamTest") << "DataSource is NOT assigned!";
        return;
    }

    // Common intialize
    DATA_DIR_PATH_ = "/home/towa/alice3/hist/"; // 自動で取得できるようにしたい
    TIME_ = plot_histogram::currentDateTime();
    canvas_  = new TCanvas("canvas", "canvas", 800, 600);
    gStyle->SetOptStat(0);

    title_latex_.SetTextSize(0.04);
    title_latex_.SetTextFont(62);
    condition_latex_.SetTextSize(0.03);
    condition_latex_.SetTextFont(62);
}

plot_BeamTest::~plot_BeamTest() {
}

// Plot for KEK202412 results
void plot_BeamTest::run_kek_plots() {
    if (source_ != DataSource::KEK202412) {
        LOG_ERROR.source("plot_BeamTest::run_kek_plots") << "This method is for KEKE data.";
        return;
    }
    LOG_STATUS.source("plot_BeamTest::run_kek_plots") << "Starting KEK plot generation.";

    // plot1--Cluster Charge (blk, std, gap)
    std::vector<PlotConfig> threshold_configs = {
        {"22p5", "std", "10", "500",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"},
         "22.5um, std, 10V", kRed, 20, 1},
        {"22p5", "std", "7", "500",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"},
         "22.5um, std, 7V", kRed+2, 20, 1},
        {"22p5", "std", "4", "500",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"},
         "22.5um, std, 4V", kRed+3, 20, 1},
        {"22p5", "blk", "10", "300",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300"},
         "22.5um, blk, 10V", kBlack, 21, 2},
        // {"22p5", "blk", "7", "300",
        //  {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300"},
        //  "22.5um, blk, 7V", kBlack, 21, 2},
        // {"22p5", "blk", "4", "300",
        //  {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300"},
        //  "22.5um, blk, 4V", kBlack, 21, 2},
        {"22p5", "gap", "10", "500",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"},
         "22.5um, gap, 10V", kBlue, 22, 3},
        {"22p5", "gap", "7", "500",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"},
         "22.5um, gap, 7V", kBlue+2, 22, 3},
        {"22p5", "gap", "4", "500",
         {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"},
         "22.5um, gap, 4V", kBlue+3, 22, 3},
    };

    std::vector<TGraphErrors*> kek_reso_graphs;
    for(const auto& conf : threshold_configs) {
        kek_reso_graphs.push_back(create_graph_data(conf, "resolution"));
    }

    draw_multigraph("kek_resolution",
                    ";threshold [ADC];resolution in x [um]",
                    Form("plot/%s/kek_resolution_thresholdScan.pdf", NAME_.c_str()),
                    threshold_configs,
                    kek_reso_graphs,
                    {3, 10},
                    {40, 610});

    std::vector<TGraphErrors*> kek_clsize_graphs;
    for(const auto& conf : threshold_configs) {
        kek_clsize_graphs.push_back(create_graph_data(conf, "clustersize"));
    }

    draw_multigraph("kek_clustersize",
                    ";threshold [ADC];mean cluster size",
                    Form("plot/%s/kek_clustersize_thresholdScan.pdf", NAME_.c_str()),
                    threshold_configs,
                    kek_clsize_graphs,
                    {1, 10},
                    {40, 610});

    LOG_STATUS.source("plot_BeamTest::run_kek_plots") << "Creating residual check plots.";

    if(!threshold_configs.empty()) {
        //check_residual_fits(threshold_configs[0], Form("plot/%s/residuals/check_residuals_std_10V.pdf", NAME_.c_str()));
        for(const auto& config_to_check : threshold_configs) {
            std::string safe_label = config_to_check.legend_label;

            std::replace(safe_label.begin(), safe_label.end(), ' ', '_');
            //sstd::replace(safe_label.begin(), safe_label.end(), ',', ""); // 空白に変換できない
            std::replace(safe_label.begin(), safe_label.end(), '.', 'p');
            safe_label.erase(std::remove(safe_label.begin(), safe_label.end(), ','), safe_label.end());

            std::string output_filename = Form("plot/%s/residuals/check_residuals_%s.pdf", NAME_.c_str(), safe_label.c_str());

            check_residual_fits(config_to_check, output_filename);
        }
    }
}

void plot_BeamTest::run_sps_plots() {
    if(source_ != DataSource::SPS202404) {
        LOG_ERROR.source("plot_BeamTest::run_sps_plots") << "This method is for SPS data.";
        return;
    }

    LOG_STATUS.source("plot_BeamTest::run_sps_plots") << "Starting SPS plot generation.";

    std::vector<PlotConfig> threshold_configs = {
        {"15", "std", "10", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "15um, std, 10V", kRed, 20, 1},
        {"15", "std", "4", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "15um, std, 4V", kRed+3, 20, 1},
        {"22p5", "std", "10", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "22.5um, std, 10V", kRed, 21, 1},
        {"22p5", "std", "4", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "22.5um, std, 4V", kRed+3, 21, 1},
        {"15", "gap", "10", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "15um, gap, 10V", kBlue, 20, 1},
        {"15", "gap", "4", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "15um, gap, 4V", kBlue+3, 20, 1},
        {"22p5", "gap", "10", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "22.5um, gap, 10V", kBlue, 21, 1},
        {"22p5", "gap", "4", "1000",
         {"200", "300", "400", "500", "600", "700", "800", "900", "1000"},
         "22.5um, gap, 4V", kBlue+3, 21, 1},
    };

    std::vector<TGraphErrors*> sps_reso_graphs;
    for(const auto& conf : threshold_configs) {
        sps_reso_graphs.push_back(create_graph_data(conf ,"resolution"));
    }
    draw_multigraph("sps_resolution",
                    ";threshold [ADC];resolution in x [um]",
                    Form("plot/%s/sps_resolution_thresholdScan.pdf", NAME_.c_str()),
                    threshold_configs,
                    sps_reso_graphs,
                    {0, 10},
                    {190, 1010});

    std::vector<TGraphErrors*> sps_clsize_graphs;
    for(const auto& conf : threshold_configs) {
        sps_clsize_graphs.push_back(create_graph_data(conf, "clustersize"));
    }

    draw_multigraph("sps_clustersize",
                    ";threshold [ADC];mean cluster size",
                    Form("plot/%s/sps_clustersize_thresholdScan.pdf", NAME_.c_str()),
                    threshold_configs,
                    sps_clsize_graphs,
                    {1, 10},
                    {40, 1010});
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

    TLegend* legend_point = new TLegend(0.55, 0.6, 0.9, 0.9);
    legend_point->SetFillStyle(0);
    legend_point->SetBorderSize(0);
    legend_point->SetTextSize(0.04);

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
    std::vector<double> x_vals, y_vals, x_errs, y_errs;

    for(const auto& thd : config.scan_values) {
        std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
                                           DATA_DIR_PATH_.c_str(),
                                           NAME_.c_str(),
                                           NAME_.c_str(),
                                           config.pixel_pitch.c_str(),
                                           config.chip_type.c_str(),
                                           config.voltage.c_str(),
                                           config.seed_thd.c_str(),
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
    canvas_->SetLeftMargin(0.13);
    canvas_->SetRightMargin(0.07);

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(title.c_str());

    TLegend* legend_point = new TLegend(0.50, 0.65, 0.90, 0.91);
    legend_point->SetFillStyle(0);
    legend_point->SetBorderSize(0);
    legend_point->SetTextSize(0.035);
    legend_point->SetNColumns(2);

    TLegend* legend_line = new TLegend(0.50, 0.65, 0.90, 0.91);
    legend_line->SetFillStyle(0);
    legend_line->SetBorderSize(0);
    legend_line->SetTextSize(0.035);
    legend_line->SetTextColor(kWhite);
    legend_line->SetNColumns(2);

    for(size_t i=0; i<configs.size(); ++i) {
        if(!graphs[i] || graphs[i]->GetN() == 0 || graphs[i] == nullptr) {
            LOG_WARNING.source("plot_BeamTest::draw_multigraph") << "Graph not found.";
            continue; 
        }

        graphs[i]->SetLineColor(configs[i].color);
        graphs[i]->SetMarkerColor(configs[i].color);
        graphs[i]->SetMarkerStyle(configs[i].marker_style);
        graphs[i]->SetLineStyle(configs[i].line_style);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerSize(1.1);

        mg->Add(graphs[i], "L");
        mg->Add(graphs[i], "PE");
        legend_line->AddEntry(graphs[i], configs[i].legend_label.c_str(), "l");
        legend_point->AddEntry(graphs[i], configs[i].legend_label.c_str(), "pe");
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
    condition_latex_.DrawLatexNDC(0.15, 0.91, BEAM_INFO_.c_str());
    condition_latex_.DrawLatexNDC(0.68, 0.91, Form("Plotted on %s", TIME_.c_str()));
    //title_latex_.DrawLatexNDC(0.15, 0.83, canvas_title.c_str());

    canvas_->SaveAs(output_filename.c_str());
}

void plot_BeamTest::check_residual_fits(
    const PlotConfig& config,
    const std::string& output_filename) {
    LOG_STATUS.source("plot_BeamTest::check_residual_fits") << "Creating residual check plots for " << config.legend_label;

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
            std::string input_file_path = Form("%s%s/%s_%s_%s_%sV_SeedThd%se_NeighborThd%se.root",
                                               DATA_DIR_PATH_.c_str(),
                                               NAME_.c_str(),
                                               NAME_.c_str(),
                                               config.pixel_pitch.c_str(),
                                               config.chip_type.c_str(),
                                               config.voltage.c_str(),
                                               config.seed_thd.c_str(),
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