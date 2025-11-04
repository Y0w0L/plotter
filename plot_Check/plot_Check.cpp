#include "plot_Check.h"

plot_Check::plot_Check() {
    LOG_STATUS.source("plot_Check::plot_Check") << "plot_Check object is created.";
}

plot_Check::~plot_Check() {
    cleanupResources();
}

TFile* plot_Check::openFile(const std::string& fileName) {
    if(m_fileCache.find(fileName) == m_fileCache.end()) {
        TFile* file = TFile::Open(fileName.c_str());
        if(!file || file->IsZombie()) {
            LOG_ERROR.source("plot_Check::openFile") << "Failed to open file: " << fileName;
            delete file;
            m_fileCache[fileName] = nullptr;
        } else {
            LOG_STATUS.source("plot_Check::openFile") << "Opend and cached file: " << fileName;
            m_fileCache[fileName] = file;
        }
    }
    return m_fileCache[fileName];
}

void plot_Check::cleanupResources() {
    for (auto const& [key, file] : m_fileCache) {
        if(file) {
            file->Close();
            delete file;
        }
    }
    m_fileCache.clear();
    LOG_STATUS.source("plot_Check::cleanupResources") << "Cleaned up file cache.";
}

void plot_Check::setHistStyle(TH1D* hist, int color, int marker, float alpha) {
    if(!hist) return;
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(marker);
    hist->SetMarkerSize(0.8);
    hist->SetFillColorAlpha(color, alpha);
}

// TH1D* plot_Check::getScaledHist(const std::string& file, const std::string& histName,
//                                 double scaleFactor, const std::string& title, bool isMerged) {
//     TH1D* hist = nullptr;
//     if(isMerged) {
//         hist = plot_BeamTest::get_merged_object<TH1D>(file, histName);
//     } else {
//         TFile* f = openFile(file);
//         if(f) {
//             TH1* h_obj = (TH1*)f->Get(histName.c_str());
//             if(h_obj) {
//                 hist = (TH1D*)h_obj->Clone();
//             }
//         }
//     }
// }

TH1D* plot_Check::scalingHistogram(TH1D* hist, double scale_factor, std::string title) {
    if (!hist) {
        LOG_WARNING.source("plot_Check::scalignHistogram") << "scalingHistogram received a null histogram." << std::endl;
        return nullptr;
    }

    LOG_STATUS << "Scaling histogram: " << hist->GetName() << " by factor " << scale_factor << std::endl;

    TAxis* xAxis = hist->GetXaxis();
    int nbins = xAxis->GetNbins();
    double xmin = xAxis->GetXmin();
    double xmax = xAxis->GetXmax();

    // 新しいX軸の範囲を計算
    double new_xmin = xmin * scale_factor;
    double new_xmax = xmax * scale_factor;

    // ヒストグラムの軸を直接設定 (in-place)
    xAxis->Set(nbins, new_xmin, new_xmax);
    
    // タイトルと軸ラベルを更新
    hist->SetTitle(title.c_str());

    // ビンの中身をチェック (元のコードから流用)
    for (int i = 0; i <= nbins + 1; ++i) { 
        double content_orig = hist->GetBinContent(i);
        if (content_orig < 0) {
            hist->SetBinContent(i, 0);
        }
    }

    return hist; // 元のヒストグラムへのポインタを返す
}

TH1D* plot_Check::getScaledHist(const std::string& file, const std::string& histName, 
                                double scaleFactor, const std::string& title, bool isMerged) 
{
    TH1D* hist = nullptr;
    if (isMerged) {
        hist = plot_BeamTest::get_merged_object<TH1D>(file, histName);
    } else {
        TFile* f = openFile(file);
        if (f) {
            TH1* h_obj = (TH1*)f->Get(histName.c_str());
            if (h_obj) {
                // クローンを作成 (重要：元のファイル内のヒストグラムを変更しないため)
                hist = (TH1D*)h_obj->Clone(Form("%s_clone", h_obj->GetName())); 
            }
        }
    }

    if (!hist) {
        LOG_WARNING.source("plot_Check::getScaledHist") << "Could not get hist: " << histName << " from file: " << file << std::endl;
        return nullptr;
    }
    
    // scalingHistogram は渡されたヒストグラム(hist)を直接変更する
    TH1D* scaledHist = this->scalingHistogram(hist, scaleFactor, title);
    
    // delete hist; // <--- この行を削除、またはコメントアウトします
    
    return scaledHist; // scaledHist と hist は同じオブジェクトを指している
}


void plot_Check::run_plotCheck() {
    LOG_STATUS.source("plot_Check::run_plotCheck") << "Start run_plotCheck";

    std::vector<DatasetConfig> datasets;

    datasets.push_back({
        "SQ P15 GAP 10V", "_gap",
        "/home/towa/alice3/hist/sps_check/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e", 0.238,
        "/home/towa/alice3/plotter/tools/analysis_py_sq_p15_gap_10v_n0e_st0_nt0_pip_120GeV_masetti_cps1.root", 1000.0
    });

    datasets.push_back({
        "SQ P15 STD 10V", "_std",
        "/home/towa/alice3/hist/sps_check/sps202404_15_std_10V_SeedThd1000e_NeighborThd200e", 0.240,
        "/home/towa/alice3/plotter/tools/analysis_py_sq_p15_std_10v_n0e_st0_nt0_pip_120GeV_masetti_cps1.root", 1000.0
    });

    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    TLegend* legend = new TLegend(0.5, 0.65, 0.9, 0.8);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);

    for(const auto& ds : datasets) {
        LOG_STATUS << "Processing dataset: " << ds.name;

        plotAllCharge(canvas, legend, ds);

        plotClusterSize(canvas, legend, ds);

        plotComparison(canvas, legend, ds, 
                       "AnalysisCE65/CE65_6/cluster/clusterSeedCharge", "seed_charge",
                       "exp, seed", "sim, seed",
                       "clusterSeedCharge_sim_exp", ds.name,
                       10, 0, 4000);
        
        plotComparison(canvas, legend, ds,
                       "AnalysisCE65/CE65_6/cluster/clusterCharge", "cluster_charge",
                       "exp, cluster", "sim, cluster",
                       "clusterCharge_sim_exp", ds.name,
                       10, 0, 4000);

        plotComparison(canvas, legend, ds,
                       "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum", "cluster_neighbor_charge_sum",
                       "exp, neighbor sum", "sim, neighbor sum",
                       "clusterNeighborChargeSum_sim_exp", ds.name,
                       10, 0, 2000);

        plotSeedChargeByCS(canvas, legend, ds);
    }

    delete canvas;
    delete legend;
    cleanupResources();

    LOG_STATUS.source("plot_Check::run_plotCheck") << "Finished run_plotCheck";
}

void plot_Check::plotAllCharge(TCanvas* c, TLegend* l, const DatasetConfig& ds) {
    LOG_STATUS.source("plot_Check::plotAllCharge") << "Plotting all charge for " << ds.name; 
    c->Clear();
    l->Clear();

    TH1D* h_exp_cl = getScaledHist(ds.expFileBase, "AnalysisCE65/CE65_6/cluster/clusterCharge", ds.expScale, ";charge [e];counts", true);
    TH1D* h_exp_sd = getScaledHist(ds.expFileBase, "AnalysisCE65/CE65_6/cluster/clusterSeedCharge", ds.expScale, ";charge [e];counts", true);
    TH1D* h_sim_cl = getScaledHist(ds.simFile, "cluster_charge", ds.simScale, ";charge [e];counts");
    TH1D* h_sim_sd = getScaledHist(ds.simFile, "seed_charge", ds.simScale, ";charge [e];counts");

    if (!h_exp_cl || !h_exp_sd || !h_sim_cl || !h_sim_sd) {
         LOG_ERROR.source("plot_Check::plotAllCharge") << "Missing histograms for plotAllCharge";
         delete h_exp_cl; delete h_exp_sd; delete h_sim_cl; delete h_sim_sd;
         return;
    }

    setHistStyle(h_exp_cl, kBlue, 20);
    setHistStyle(h_exp_sd, kAzure - 3, 24);
    setHistStyle(h_sim_cl, kRed + 1, 20);
    setHistStyle(h_sim_sd, kPink - 2, 24);

    h_exp_cl->Rebin(10);
    h_exp_sd->Rebin(10);

    if(h_exp_cl->GetMaximum() > 0) h_exp_cl->Scale(1.0 / h_exp_cl->GetMaximum());
    if(h_exp_sd->GetMaximum() > 0) h_exp_sd->Scale(1.0 / h_exp_sd->GetMaximum());
    if(h_sim_cl->GetMaximum() > 0) h_sim_cl->Scale(1.0 / h_sim_cl->GetMaximum());
    if(h_sim_sd->GetMaximum() > 0) h_sim_sd->Scale(1.0 / h_sim_sd->GetMaximum());

    h_sim_cl->GetXaxis()->SetRangeUser(0, 4000);
    h_sim_cl->GetYaxis()->SetRangeUser(0, 1.1);

    h_sim_cl->Draw("PE");
    h_exp_cl->Draw("samePE");
    h_sim_sd->Draw("samePE");
    h_exp_sd->Draw("samePE");

    l->AddEntry(h_exp_sd, "exp, seed", "pe");
    l->AddEntry(h_exp_cl, "exp, cluster", "pe");
    l->AddEntry(h_sim_sd, "sim, seed", "pe");
    l->AddEntry(h_sim_cl, "sim, cluster", "pe");
    l->Draw();

    TLatex title;
    title.SetTextAlign(12);
    title.SetTextSize(0.05);
    title.DrawLatexNDC(0.6, 0.85, ds.name.c_str());

    c->SaveAs(Form("./plot/AllclusterCharge_sim_exp%s.pdf", ds.suffix.c_str()));
    
    // 5. Clean up
    delete h_exp_cl; delete h_exp_sd; delete h_sim_cl; delete h_sim_sd;
}

void plot_Check::plotClusterSize(TCanvas* c, TLegend* l,  const DatasetConfig& ds) {
    LOG_STATUS.source("plot_Check::plotClusterSize") << "Plotting ClusterSize for " << ds.name;
    c->Clear();
    l->Clear();

    TH1D* h_exp_size = plot_BeamTest::get_merged_object<TH1D>(ds.expFileBase, "AnalysisCE65/CE65_6/cluster/clusterSize");
    
    TH1D* h_sim_size_orig = nullptr;
    TFile* f_sim = openFile(ds.simFile);
    if(f_sim) h_sim_size_orig = (TH1D*)f_sim->Get("cluster_size");

    if (!h_exp_size || !h_sim_size_orig) {
         LOG_ERROR.source("plot_Check::plotClsuterSize") << "Missing histograms for plotClusterSize";
         delete h_exp_size;
         return;
    }
    h_exp_size = (TH1D*)h_exp_size->Clone();
    TH1D* h_sim_size = (TH1D*)h_sim_size_orig->Clone();

    setHistStyle(h_exp_size, kBlue, 20, 0.2f);
    setHistStyle(h_sim_size, kRed + 1, 20, 0.2f);

    if(h_exp_size->GetEntries() > 0) h_exp_size->Scale(1.0 / h_exp_size->GetEntries());
    if(h_sim_size->GetEntries() > 0) h_sim_size->Scale(1.0 / h_sim_size->GetEntries());
    
    h_sim_size->GetXaxis()->SetRangeUser(1, 10);
    h_sim_size->GetYaxis()->SetRangeUser(0, 1);

    h_sim_size->Draw("HIST");
    h_exp_size->Draw("same HIST");
    h_sim_size->Draw("same PE");
    h_exp_size->Draw("samePE");

    l->AddEntry(h_exp_size, "exp", "pef");
    l->AddEntry(h_exp_size, Form("mean = %f", h_exp_size->GetMean()), "");
    l->AddEntry(h_sim_size, "sim", "pef");
    l->AddEntry(h_sim_size, Form("mean = %f", h_sim_size->GetMean()), "");
    l->Draw();

    TLatex title;
    title.SetTextAlign(12);
    title.SetTextSize(0.05);
    title.DrawLatexNDC(0.6, 0.85, ds.name.c_str());
    
    c->SaveAs(Form("./plot/clusterSize_sim_exp%s.pdf", ds.suffix.c_str()));

    delete h_exp_size;
    delete h_sim_size;
}

void plot_Check::plotComparison(TCanvas* c, TLegend* l, const DatasetConfig& ds, 
                                const std::string& expHistName, const std::string& simHistName,
                                const std::string& expLegend, const std::string& simLegend,
                                const std::string& outName, const std::string& plotTitle,
                                int rebin, double xMin, double xMax, bool normalizeToMax) {
    LOG_STATUS.source("plot_Check::plotComparison") << "Plotting Comparison " << outName << " for " << ds.name;
    c->Clear();
    l->Clear();

    std::string title = ";charge [e];counts";
    TH1D* h_exp = getScaledHist(ds.expFileBase, expHistName, ds.expScale, title, true);
    TH1D* h_sim = getScaledHist(ds.simFile, simHistName, ds.simScale, title);

    if (!h_exp || !h_sim) {
         LOG_ERROR << "Missing histograms for " << outName << std::endl;
         delete h_exp; delete h_sim;
         return;
    }

    setHistStyle(h_exp, kBlue, 20);
    setHistStyle(h_sim, kRed + 1, 20);

    if (rebin > 1) {
        h_exp->Rebin(rebin);
    }

    if (normalizeToMax) {
        if(h_exp->GetMaximum() > 0) h_exp->Scale(1.0 / h_exp->GetMaximum());
        if(h_sim->GetMaximum() > 0) h_sim->Scale(1.0 / h_sim->GetMaximum());
    }

    h_sim->GetXaxis()->SetRangeUser(xMin, xMax);

    h_sim->Draw("PE");
    h_exp->Draw("samePE");

    l->AddEntry(h_exp, expLegend.c_str(), "pe");
    l->AddEntry(h_sim, simLegend.c_str(), "pe");
    l->Draw();
    
    TLatex latexTitle;
    latexTitle.SetTextAlign(12);
    latexTitle.SetTextSize(0.05);
    latexTitle.DrawLatexNDC(0.6, 0.85, plotTitle.c_str());

    c->SaveAs(Form("./plot/%s%s.pdf", outName.c_str(), ds.suffix.c_str()));
    
    delete h_exp;
    delete h_sim;
}

void plot_Check::plotSeedChargeByCS(TCanvas* c, TLegend* l, const DatasetConfig& ds) {
    LOG_STATUS.source("plot_Check::plotSeedChargeByCS") << "Plotting SeedChargeByCS for " << ds.name;

    for(int cs = 1; cs <= 6; ++cs) {
        c->Clear();
        l->Clear();

        std::string expHistName = Form("AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size%d", cs);
        std::string simHistName = Form("seed_charge_size_%d", cs);
        std::string plotTitle = Form("%s cs%d", ds.name.c_str(), cs);
        std::string outName = Form("seed_charge_cs/clusterSeedCharge_cs%d_sim_exp", cs);

        std::string title = ";charge [e];counts";
        TH1D* h_exp = getScaledHist(ds.expFileBase, expHistName, ds.expScale, title, true);
        TH1D* h_sim = getScaledHist(ds.simFile, simHistName, ds.simScale, title);

        if (!h_exp || !h_sim) {
            LOG_ERROR << "Missing histograms for " << outName << std::endl;
            delete h_exp; delete h_sim;
            continue; // next cluster size
        }

        setHistStyle(h_exp, kBlue, 24);
        setHistStyle(h_sim, kRed, 24);
        
        h_exp->Rebin(10);

        if (h_sim->GetMaximum() > 0) h_sim->Scale(1.0 / h_sim->GetMaximum());
        if (h_exp->GetMaximum() > 0) h_exp->Scale(1.0 / h_exp->GetMaximum());

        h_sim->GetXaxis()->SetRangeUser(0, 4000);
        
        h_sim->Draw("PE");
        h_exp->Draw("samePE");

        l->AddEntry(h_exp, "exp, seed", "pe");
        l->AddEntry(h_sim, "sim, seed", "pe");
        l->Draw();

        TLatex latexTitle;
        latexTitle.SetTextAlign(12);
        latexTitle.SetTextSize(0.05);
        latexTitle.DrawLatexNDC(0.5, 0.85, plotTitle.c_str());

        c->SaveAs(Form("./plot/%s%s.pdf", outName.c_str(), ds.suffix.c_str()));
        
        delete h_exp;
        delete h_sim;
    }
}

void plot_Check::plotDepositedCharge(TCanvas* c, TLegend* l) {
    LOG_STATUS.source("plot_Check::plotDepositedCharge") << "Plotting DepositedCharge";
    c->Clear();
    l->Clear();

    std::string file = "/home/towa/alice3/hist/test/CE65_sq_p15_gap_10v_pip_120GeV_default.root";
    
    TH1D* h_dep = getScaledHist(file, "DepositionGeant4/deposited_charge_CE65", 1000, ";charge [e];counts");
    TH1D* h_cls = getScaledHist(file, "DetectorHistogrammer/CE65/charge/cluster_charge", 1000, ";charge [e];counts");
    
    if (!h_dep || !h_cls) {
         LOG_ERROR << "Missing histograms for plotDepositedCharge" << std::endl;
         delete h_dep; delete h_cls;
         return;
    }

    setHistStyle(h_dep, kBlue, 20);
    setHistStyle(h_cls, kRed + 1, 20);
    
    h_dep->GetXaxis()->SetRangeUser(0, 20000);
    
    h_dep->Draw("PE");
    h_cls->Draw("samePE");

    l->AddEntry(h_dep, "deposited charge", "pe");
    l->AddEntry(h_dep, Form("mean = %f", h_dep->GetMean()), "");
    l->AddEntry(h_cls, "cluster charge", "pe");
    l->AddEntry(h_cls, Form("mean = %f", h_cls->GetMean()), "");
    l->Draw();

    c->SaveAs("./plot/depositedCharge.pdf");
    
    delete h_dep;
    delete h_cls;
}

// void plot_Check::run_Check() {
//     LOG_STATUS.source("plot_Check::run_Check") << "Start run_Check";

//     //plot_BeamTest plot_BeamTest;

//     // std::string input_filename_sim = Form("/home/towa/alice3/hist/ce65sim202505/ce65sim202505_15_gap_10V_SeedThd0e_NeighborThd25e_0.root");
//     // std::string input_filename_exp = Form("/home/towa/alice3/hist/sps202404/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e_0.root");

//     //std::string input_filename_sim = Form("/home/towa/alice3/plotter/tools/analysis_py_p15_gap_sq_10v_st0_nt60_n60e_2um.root");
//     //std::string input_filename_exp = Form("/home/towa/alice3/hist/sps202404/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e_0.root");

//     //TFile* inputROOTFile_sim = TFile::Open(input_filename_sim.c_str());
//     //TFile* inputROOTFile_exp = TFile::Open(input_filename_exp.c_str());

//     //TH1D* h_clusterCharge_sim = (TH1D*)inputROOTFile_sim->Get("seed_charge");
//     //TH1D* h_clusterCharge_exp = (TH1D*)inputROOTFile_exp->Get("AnalysisCE65/CE65_6/cluster/clusterCharge");

//     std::string hist_name = ("AnalysisCE65/CE65_6/cluster/clusterSeedCharge");

//     std::string base_file_path = "/home/towa/alice3/hist/sps202404/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e";
//     TH1D* h_clusterCharge_exp = plot_BeamTest::get_merged_object<TH1D>(base_file_path, hist_name);

//     double scale_factor = 0.238; // SQ P15 GAP 10V
//     //double scale_factor = 0.240; // SQ P15 STD 10V
//     //double scale_factor = 0.215; // SQ P22.5 GAP 10V
//    // double scale_factor = 0.235; // SQ P22.5 STD 10V
//    //double scale_factor = 1;

//     TAxis* xAxis_orig = h_clusterCharge_exp->GetXaxis();
//     int nbins = xAxis_orig->GetNbins();
//     double xmin = xAxis_orig->GetXmin();
//     double xmax = xAxis_orig->GetXmax();

//     double new_xmin = xmin * scale_factor;
//     double new_xmax = xmax * scale_factor;

//     TH1D* h_clusterCharge_exp_scaled = new TH1D("h_clusterCharge_exp_scaled", "", nbins, 0, new_xmax);
//     for(int i=0; i <= nbins + 1; ++i) {
//         double content = h_clusterCharge_exp->GetBinContent(i) * scale_factor;
//         double error = h_clusterCharge_exp->GetBinError(i) * scale_factor;
//         h_clusterCharge_exp_scaled->SetBinContent(i, content);
//         h_clusterCharge_exp_scaled->SetBinError(i, error);
//     }

// //    double sim_x_scale_factor = 1000.0; // X軸の値を1000倍する
// double sim_x_scale_factor = 1000;

//     //TAxis* xAxis_sim_orig = h_clusterCharge_sim->GetXaxis();
//     // int sim_nbins = xAxis_sim_orig->GetNbins();
//     // double sim_xmin = xAxis_sim_orig->GetXmin();
//     // double sim_xmax = xAxis_sim_orig->GetXmax();

//     // 新しいX軸の範囲を計算（ke -> eなので、最小値と最大値を1000倍）
//     // double sim_new_xmin = sim_xmin * sim_x_scale_factor;
//     // double sim_new_xmax = sim_xmax * sim_x_scale_factor;

//     // 新しいX軸を持つヒストグラムを新規作成
//     // TH1D* h_clusterCharge_sim_scaled = new TH1D("h_clusterCharge_sim_scaled",
//     //                                             "Cluster Seed Charge;charge [e];counts", // X軸ラベルを "charge [e]" に変更
//     //                                             sim_nbins, sim_new_xmin, sim_new_xmax);


//     // for (int i = 0; i <= sim_nbins + 1; ++i) { // アンダーフロー/オーバーフロービンも考慮
//     //     double content_orig = h_clusterCharge_sim->GetBinContent(i);
//     //     double error_orig   = h_clusterCharge_sim->GetBinError(i);

//     //     // 負のビンがある場合、0に置き換える
//     //     if (content_orig < 0) {
//     //         content_orig = 0;
//     //     }

//     //     h_clusterCharge_sim_scaled->SetBinContent(i, content_orig);
//     //     h_clusterCharge_sim_scaled->SetBinError(i, error_orig);
//     // }

//     //double sim_y_scale_factor = 1000.0; // Y軸を1000倍する目的
//     // double sim_y_scale_factor = 1;

//     // for (int i = 0; i <= sim_nbins + 1; ++i) {
//     //     h_clusterCharge_sim_scaled->SetBinContent(i, h_clusterCharge_sim_scaled->GetBinContent(i) * sim_y_scale_factor);
//     //     h_clusterCharge_sim_scaled->SetBinError(i, h_clusterCharge_sim_scaled->GetBinError(i) * sim_y_scale_factor);
//     // }

//     // h_clusterCharge_sim_scaled->SetLineColor(kRed);
//     // //h_clusterCharge_exp_scaled->Rebin(10);
//     // h_clusterCharge_exp_scaled->SetLineColor(kBlue);

//     // h_clusterCharge_sim_scaled->SetMarkerSize(0.8);
//     // h_clusterCharge_sim_scaled->SetMarkerStyle(20);
//     // h_clusterCharge_sim_scaled->SetMarkerColor(kRed);

//     // h_clusterCharge_exp_scaled->SetMarkerSize(0.8);
//     // h_clusterCharge_exp_scaled->SetMarkerStyle(20);
//     // h_clusterCharge_exp_scaled->SetMarkerColor(kBlue);

//     // h_clusterCharge_sim_scaled->Scale(1.0/h_clusterCharge_sim_scaled->GetMaximum());
//     // h_clusterCharge_exp_scaled->Scale(1.0/h_clusterCharge_exp_scaled->GetMaximum());
//     // h_clusterCharge_sim_scaled->Scale(1.0/h_clusterCharge_sim->GetEntries());
//     // h_clusterCharge_exp_scaled->Scale(1.0/h_clusterCharge_exp->GetEntries());

//     // TCanvas* c = new TCanvas("c", "c", 800, 600);

//     // TLegend* legend = new TLegend(0.5, 0.65, 0.9, 0.8);
//     // legend->SetFillStyle(0);
//     // legend->SetBorderSize(0);
//     // legend->SetTextSize(0.03);

//     // TLatex title;
//     // title.SetTextAlign(12);
//     // title.SetTextSize(0.05);
//     // //title.DrawLatexNDC(0.1, 0.85, "SQ P15 GAP 10V");

//     // legend->AddEntry(h_clusterCharge_sim_scaled, "Simulation", "pe");
//     // legend->AddEntry(h_clusterCharge_exp_scaled, "Experiment", "pe");

//     // h_clusterCharge_sim_scaled->GetXaxis()->SetRangeUser(0, 4000);
//     // h_clusterCharge_sim_scaled->GetYaxis()->SetRangeUser(0, 1.1);

//     // h_clusterCharge_sim_scaled->Draw("histPE");
//     // h_clusterCharge_exp_scaled->Draw("samehistPE");
//     // title.DrawLatexNDC(0.6, 0.85, "SQ P15 GAP 10V");
//     // legend->Draw();

//     // c->SaveAs("./plot/check.pdf");
// }

// TH1D* plot_Check::scalingHistogram(TH1D* hist, double scale_factor, std::string title) {
//     LOG_STATUS.source("plot_Check::scalingHistogram") << "Start scalingHistogram";

//     TAxis* xAxis = hist->GetXaxis();
//     int nbins = xAxis->GetNbins();
//     double xmin = xAxis->GetXmin();
//     double xmax = xAxis->GetXmax();

//     // 新しいX軸の範囲を計算（ke -> eなので、最小値と最大値を1000倍）
//     double new_xmin = xmin * scale_factor;
//     double new_xmax = xmax * scale_factor;

//     // 新しいX軸を持つヒストグラムを新規作成
//     TH1D* scaled_hist = new TH1D("hist",
//                           title.c_str(),
//                           nbins,
//                           new_xmin,
//                           new_xmax);

//     double content_orig, error_orig;

//     for (int i = 0; i <= nbins + 1; ++i) { // アンダーフロー/オーバーフロービンも考慮
//         content_orig = hist->GetBinContent(i);
//         error_orig   = hist->GetBinError(i);

//         // 負のビンがある場合、0に置き換える
//         if (content_orig < 0) {
//             content_orig = 0;
//         }

//         scaled_hist->SetBinContent(i, content_orig);
//         scaled_hist->SetBinError(i, error_orig);
//     }

//     return scaled_hist;
// }

// TH1D* plot_Check::scalingClusterSize(TH1D* hist) {
//     int nbins = hist->GetNbinsX();
//     double xlow = 0.5;
//     double xup = 0.5 + nbins;
//     TAxis* xaxis = hist->GetXaxis();
//     xaxis->Set(nbins, xlow, xup);
//     return hist;
// }

// void plot_Check::run_plotCheck() {
//     LOG_STATUS.source("plot_Check::run_plotCheck") << "Start run_plotCheck";

//     // Scale factor
//     double scale_factor_sq_p15_gap_10v = 0.238; // SQ P15 GAP 10V
//     double scale_factor_sq_p15_std_10v = 0.240; // SQ P15 STD 10V
//     double scale_factor_sq_p225_gap_10v = 0.215; // SQ P22.5 GAP 10V
//     double scale_factor_sq_p225_std_10v = 0.235; // SQ P22.5 STD 10V
//     double scale_factor_charge_sim = 1000;

//     // filename
//     std::string filename_exp_sq_p15_gap_10v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e";
//     std::string filename_exp_sq_p15_std_10v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_15_std_10V_SeedThd1000e_NeighborThd200e";
//     std::string filename_exp_sq_p225_gap_10v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_22p5_gap_10V_SeedThd1000e_NeighborThd200e";
//     std::string filename_exp_sq_p225_std_10v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_22p5_std_10V_SeedThd1000e_NeighborThd200e";

//     std::string filename_exp_sq_p15_gap_4v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_15_gap_4V_SeedThd1000e_NeighborThd200e";
//     std::string filename_exp_sq_p15_std_4v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_15_std_4V_SeedThd1000e_NeighborThd200e";
//     std::string filename_exp_sq_p225_gap_4v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_22p5_gap_4V_SeedThd1000e_NeighborThd200e";
//     std::string filename_exp_sq_p225_std_4v_st1000_nt200 = "/home/towa/alice3/hist/sps_check/sps202404_22p5_std_4V_SeedThd1000e_NeighborThd200e";

//     //std::string filename_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = "/home/towa/alice3/plotter/tools/analysis_py_p15_gap_sq_10v_st0_nt60.root";
//     // std::string filename_sim_sq_p15_gap_10v_st0_nt60_n60e_2um = "/home/towa/alice3/plotter/tools/analysis_py_p15_gap_sq_10v_st0_nt60_n60e_2um.root";
//     // std::string filename_sim_sq_p15_gap_10v_st0_nt60_n60e_3um = "/home/towa/alice3/plotter/tools/analysis_py.root";
//     // std::string filename_sim_sq_p15_gap_10v_st0_nt60_n60e_4um = "";
//     // std::string filename_sim_sq_p15_gap_10v_st0_nt60_n60e_5um = "";

//     std::string filename_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = "/home/towa/alice3/plotter/tools/analysis_py_sq_p15_gap_10v_n0e_pip_120GeV_masetti_default.root";
//     std::string filename_sim_sq_p15_std_10v_st0_nt0_n0e_1um = "/home/towa/alice3/plotter/tools/analysis_py_sq_p15_std_10v_n0e_pip_120GeV_masetti_default.root";
//     std::string filename_sim_sq_p15_gap_10v_st0_nt0_n0e_1um_cps5 = "/home/towa/alice3/plotter/tools/analysis_py_sq_p15_gap_10v_n0e_pip_120GeV_masetti_cps5.root";
//     std::string filename_sim_sq_p15_std_10v_st0_nt0_n0e_1um_cps5 = "/home/towa/alice3/plotter/tools/analysis_py_sq_p15_std_10v_n0e_pip_120GeV_masetti_cps5.root";

//     // histname
//     std::string histname_clusterCharge = "AnalysisCE65/CE65_6/cluster/clusterCharge";
//     std::string histname_clusterSeedCharge = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge";
//     std::string histname_clusterSize = "AnalysisCE65/CE65_6/cluster/clusterSize";
//     std::string histname_clusterNeighborsCharge = "AnalysisCE65/CE65_6/cluster/clusterNeighborsCharge";
//     std::string histname_clusterNeighborsChargeSum = "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum";
//     std::string histname_clusterSeedCharge_clsize1 = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size1";
//     std::string histname_clusterSeedCharge_clsize2 = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size2";
//     std::string histname_clusterSeedCharge_clsize3 = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size3";
//     std::string histname_clusterSeedCharge_clsize4 = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size4";
//     std::string histname_clusterSeedCharge_clsize5 = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size5";
//     std::string histname_clusterSeedCharge_clsize6 = "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size6";

//     std::string histname_sim_clusterCharge = "cluster_charge";
//     std::string histname_sim_clusterSeedCharge = "seed_charge";
//     std::string histname_sim_clusterSize = "cluster_size";
//     std::string histname_sim_neighborChargeSum = "cluster_neighbor_charge_sum";
//     std::string histname_sim_neighborCharge = "cluster_neighbor_charge";
//     std::string histname_sim_clusterSeedCharge_size1 = "seed_charge_size_1";
//     std::string histname_sim_clusterSeedCharge_size2 = "seed_charge_size_2";
//     std::string histname_sim_clusterSeedCharge_size3 = "seed_charge_size_3";
//     std::string histname_sim_clusterSeedCharge_size4 = "seed_charge_size_4";
//     std::string histname_sim_clusterSeedCharge_size5 = "seed_charge_size_5";
//     std::string histname_sim_clusterSeedCharge_size6 = "seed_charge_size_6";

//     TH1D* h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterCharge);
//     TH1D* h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge);
//     TH1D* h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSize);
//     TH1D* h_clusterNeighborsCharge_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterNeighborsCharge);
//     TH1D* h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterNeighborsChargeSum);
//     TH1D* h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge_clsize1);
//     TH1D* h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge_clsize2);
//     TH1D* h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge_clsize3);
//     TH1D* h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge_clsize4);
//     TH1D* h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge_clsize5);
//     TH1D* h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_gap_10v_st1000_nt200, histname_clusterSeedCharge_clsize6);

//     TH1D* h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterCharge);
//     TH1D* h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge);
//     TH1D* h_clusterSize_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSize);
//     TH1D* h_clusterNeighborsCharge_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterNeighborsCharge);
//     TH1D* h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterNeighborsChargeSum);
//     TH1D* h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge_clsize1);
//     TH1D* h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge_clsize2);
//     TH1D* h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge_clsize3);
//     TH1D* h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge_clsize4);
//     TH1D* h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge_clsize5);
//     TH1D* h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200 = plot_BeamTest::get_merged_object<TH1D>(filename_exp_sq_p15_std_10v_st1000_nt200, histname_clusterSeedCharge_clsize6);

//     // TFile* inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = TFile::Open(filename_sim_sq_p15_gap_10v_st0_nt60_n60e_1um.c_str());
//     // TFile* inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_2um = TFile::Open(filename_sim_sq_p15_gap_10v_st0_nt60_n60e_1um.c_str());
//     // TFile* inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_3um = TFile::Open(filename_sim_sq_p15_gap_10v_st0_nt60_n60e_1um.c_str());

//     TFile* inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = TFile::Open(filename_sim_sq_p15_gap_10v_st0_nt0_n0e_1um.c_str());
//     TFile* inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um = TFile::Open(filename_sim_sq_p15_std_10v_st0_nt0_n0e_1um.c_str());
    
//     //TH1D* h_clusterCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_1um->Get(histname_sim_clusterCharge.c_str());
//     // TH1D* h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_1um->Get(histname_sim_clusterSeedCharge.c_str());
//     // TH1D* h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_1um->Get(histname_sim_clusterSize.c_str());

//     // TH1D* h_clusterCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_2um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_2um->Get(histname_sim_clusterCharge.c_str());
//     // TH1D* h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_2um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_2um->Get(histname_sim_clusterSeedCharge.c_str());
//     // TH1D* h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_2um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_2um->Get(histname_sim_clusterSize.c_str());

//     // TH1D* h_clusterCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_3um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_3um->Get(histname_sim_clusterCharge.c_str());
//     // TH1D* h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_3um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_3um->Get(histname_sim_clusterSeedCharge.c_str());
//     // TH1D* h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_3um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt60_n60e_3um->Get(histname_sim_clusterSize.c_str());

//     TH1D* h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterCharge.c_str());
//     TH1D* h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge.c_str());
//     TH1D* h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSize.c_str());
//     TH1D* h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_neighborChargeSum.c_str());
//     TH1D* h_clusterNeighborCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_neighborCharge.c_str());

//     TH1D* h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterCharge.c_str());
//     TH1D* h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge.c_str());
//     TH1D* h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSize.c_str());
//     TH1D* h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_neighborChargeSum.c_str());
//     TH1D* h_clusterNeighborCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_neighborCharge.c_str());

//     TH1D* h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size1.c_str());
//     TH1D* h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size2.c_str());
//     TH1D* h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size3.c_str());
//     TH1D* h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size4.c_str());
//     TH1D* h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size5.c_str());
//     TH1D* h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size6.c_str());

//     TH1D* h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size1.c_str());
//     TH1D* h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size2.c_str());
//     TH1D* h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size3.c_str());
//     TH1D* h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size4.c_str());
//     TH1D* h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size5.c_str());
//     TH1D* h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um = (TH1D*)inputROOTFile_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Get(histname_sim_clusterSeedCharge_size6.c_str());


//     // scaling
//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Cluster Charge;charge [e];counts");
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Cluster Seed Charge;charge [e];counts");
//     h_clusterNeighborsCharge_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterNeighborsCharge_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Cluster Neighbors Charge;charge [e];counts");
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Cluster Neighbors Charge Sum;charge [e];counts");

//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Cluster Charge;charge [e];counts");
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Cluster Seed Charge;charge [e];counts");
//     h_clusterNeighborsCharge_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterNeighborsCharge_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Cluster Neighbors Charge;charge [e];counts");
//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Cluster Neighbors Charge Sum;charge [e];counts");

//     // h_clusterCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = plot_Check::scalingHistogram(h_clusterCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_1um, scale_factor_charge_sim, "Cluster Charge;charge [e];counts");
//     // h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_1um  = plot_Check::scalingHistogram(h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt60_n60e_1um, scale_factor_charge_sim, "Cluster Seed Charge;charge [e];counts");

//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Cluster Charge;charge [e];counts");
//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um  = plot_Check::scalingHistogram(h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Cluster Seed Charge;charge [e];counts");
//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Cluster Charge;charge [e];counts");
//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um  = plot_Check::scalingHistogram(h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Cluster Seed Charge;charge [e];counts");
    
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Neighbor Charge Sum;charge [e];counts");
//     h_clusterNeighborCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um  = plot_Check::scalingHistogram(h_clusterNeighborCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Neighbor Charge;charge [e];counts");
//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Neighbor Charge Sum;charge [e];counts");
//     h_clusterNeighborCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um  = plot_Check::scalingHistogram(h_clusterNeighborCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Neighbor Charge;charge [e];counts");
    

//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Seed Charge clsize1;charge [e];counts");
//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Seed Charge clsize2;charge [e];counts");
//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Seed Charge clsize3;charge [e];counts");
//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Seed Charge clsize4;charge [e];counts");
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Seed Charge clsize5;charge [e];counts");
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200, scale_factor_sq_p15_gap_10v, "Seed Charge clsize6;charge [e];counts");

//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Seed Charge clsize1;charge [e];counts");
//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Seed Charge clsize2;charge [e];counts");
//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Seed Charge clsize3;charge [e];counts");
//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Seed Charge clsize4;charge [e];counts");
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Seed Charge clsize5;charge [e];counts");
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200 = plot_Check::scalingHistogram(h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200, scale_factor_sq_p15_std_10v, "Seed Charge clsize6;charge [e];counts");

//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize1;charge [e];counts");
//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize2;charge [e];counts");
//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize3;charge [e];counts");
//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize4;charge [e];counts");
//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize5;charge [e];counts");
//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize6;charge [e];counts");

//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize1;charge [e];counts");
//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize2;charge [e];counts");
//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize3;charge [e];counts");
//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize4;charge [e];counts");
//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize5;charge [e];counts");
//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingHistogram(h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um, scale_factor_charge_sim, "Seed Charge clsize6;charge [e];counts");

//     // int nbins = h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_1um->GetNbinsX();
//     // double xlow = 0.5;
//     // double xup = 0.5 + nbins;
//     // TAxis* xaxis = h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_1um->GetXaxis();
//     // xaxis->Set(nbins, xlow, xup);
//     // // cluster size
//     // h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_1um = plot_Check::scalingClusterSize(h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_1um);
//     // h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_2um = plot_Check::scalingClusterSize(h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_2um);
//     // h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_3um = plot_Check::scalingClusterSize(h_clusterSize_sim_sq_p15_gap_10v_st0_nt60_n60e_3um);

//     // h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um = plot_Check::scalingClusterSize(h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um);
//     // h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um = plot_Check::scalingClusterSize(h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um);

//     TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
//     TLegend* legend = new TLegend(0.5, 0.65, 0.9, 0.8);
//     legend->SetFillStyle(0);
//     legend->SetBorderSize(0);
//     legend->SetTextSize(0.03);

//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(20);
//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed+1);
//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed+1);
//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(20);
//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kAzure-3);
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kAzure-3);
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kPink-2);
//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kPink-2);
//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);
//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);

//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());
//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());
//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());

//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");
//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("samePE");
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200, "exp, cluster", "pe");
//     legend->AddEntry(h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->AddEntry(h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, cluster", "pe");
//     legend->Draw();

//     TLatex title;
//     title.SetTextAlign(12);
//     title.SetTextSize(0.05);
//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 GAP 10V");

//     canvas->SaveAs("./plot/AllclusterCharge_sim_exp_gap.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(20);
//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed+1);
//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed+1);
//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(20);
//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kAzure-3);
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kAzure-3);
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kPink-2);
//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kPink-2);
//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);
//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);

//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());
//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());
//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());

//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");
//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("samePE");
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200, "exp, cluster", "pe");
//     legend->AddEntry(h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->AddEntry(h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, cluster", "pe");
//     legend->Draw();

//     // TLatex title;
//     // title.SetTextAlign(12);
//     // title.SetTextSize(0.05);
//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 STD 10V");

//     canvas->SaveAs("./plot/AllclusterCharge_sim_exp_std.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(20);
//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed+1);
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed+1);
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(20);
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->GetEntries());
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetEntries());

//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(1, 10);
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetYaxis()->SetRangeUser(0, 1);

//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetFillColorAlpha(kRed-4, 0.2);
//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->SetFillColorAlpha(kBlue-4, 0.2);
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("HIST");
//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->Draw("same HIST");
//     h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("same PE");
//     h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200, "exp", "pef");
//     legend->AddEntry(h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200, Form("mean = %f", h_clusterSize_exp_sq_p15_gap_10v_st1000_nt200->GetMean()), "");
//     legend->AddEntry(h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim", "pef");
//     legend->AddEntry(h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, Form("mean = %f", h_clusterSize_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMean()), "");
//     legend->Draw();

//     // TLatex title;
//     // title.SetTextAlign(12);
//     // title.SetTextSize(0.05);
//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 GAP 10V");

//     canvas->SaveAs("./plot/clusterSize_sim_exp_gap.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(20);
//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed+1);
//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed+1);
//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(20);
//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->GetEntries());
//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetEntries());

//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(1, 10);
//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetYaxis()->SetRangeUser(0, 1);

//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetFillColorAlpha(kRed-4, 0.2);
//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->SetFillColorAlpha(kBlue-4, 0.2);

//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("HIST");
//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->Draw("same HIST");
//     h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("same PE");
//     h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->Draw("same PE");

//     legend->AddEntry(h_clusterSize_exp_sq_p15_std_10v_st1000_nt200, "exp", "fpe");
//     legend->AddEntry(h_clusterSize_exp_sq_p15_std_10v_st1000_nt200, Form("mean = %f", h_clusterSize_exp_sq_p15_std_10v_st1000_nt200->GetMean()), "");
//     legend->AddEntry(h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim", "fpe");
//     legend->AddEntry(h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um, Form("mean = %f", h_clusterSize_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMean()), "");
//     legend->Draw();

//     // TLatex title;
//     // title.SetTextAlign(12);
//     // title.SetTextSize(0.05);
//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 STD 10V");

//     canvas->SaveAs("./plot/clusterSize_sim_exp_std.pdf");

//     canvas->Clear();
//     legend->Clear();

//     TFile* inputROOTFile_sim_sq_p15_gap_10v_default = TFile::Open("/home/towa/alice3/hist/test/CE65_sq_p15_gap_10v_pip_120GeV_default.root");

//     TH1D* h_depositedCharge_CE65 = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_default->Get("DepositionGeant4/deposited_charge_CE65");
//     TH1D* h_depositedEnergy_CE65 = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_default->Get("DepositionGeant4/deposited_energy_CE65");
//     TH1D* h_clusterCharge_CE65 = (TH1D*)inputROOTFile_sim_sq_p15_gap_10v_default->Get("DetectorHistogrammer/CE65/charge/cluster_charge");

//     h_depositedCharge_CE65 = plot_Check::scalingHistogram(h_depositedCharge_CE65, 1000, ";charge [e];counts");
//     h_clusterCharge_CE65 = plot_Check::scalingHistogram(h_clusterCharge_CE65, 1000, ";charge [e];counts");

//     h_depositedCharge_CE65->SetMarkerColor(kBlue);
//     h_depositedCharge_CE65->SetLineColor(kBlue);
//     h_depositedCharge_CE65->SetMarkerStyle(20);
//     h_depositedCharge_CE65->SetMarkerSize(0.8);

//     h_clusterCharge_CE65->SetMarkerColor(kRed+1);
//     h_clusterCharge_CE65->SetLineColor(kRed+1);
//     h_clusterCharge_CE65->SetMarkerStyle(20);
//     h_clusterCharge_CE65->SetMarkerSize(0.8);

//     h_depositedCharge_CE65->GetXaxis()->SetRangeUser(0, 20000);

//     h_depositedCharge_CE65->Draw("PE");
//     h_clusterCharge_CE65->Draw("samePE");

//     legend->AddEntry(h_depositedCharge_CE65, "deposited charge", "pe");
//     legend->AddEntry(h_depositedCharge_CE65, Form("mean = %f", h_depositedCharge_CE65->GetMean()), "");
//     legend->AddEntry(h_clusterCharge_CE65, "cluster charge", "pe");
//     legend->AddEntry(h_clusterCharge_CE65, Form("mean = %f", h_clusterCharge_CE65->GetMean()), "");

//     legend->Draw();

//     canvas->SaveAs("./plot/depositedCharge.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 GAP 10V");

//     canvas->SaveAs("./plot/clusterSeedCharge_sim_exp_gap.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 STD 10V");

//     canvas->SaveAs("./plot/clusterSeedCharge_sim_exp_std.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterCharge_exp_sq_p15_gap_10v_st1000_nt200, "exp, cluster", "pe");
//     legend->AddEntry(h_clusterCharge_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, cluster", "pe");
//     legend->Draw();

//     // TLatex title;
//     // title.SetTextAlign(12);
//     // title.SetTextSize(0.05);
//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 GAP 10V");

//     canvas->SaveAs("./plot/clusterCharge_sim_exp_gap.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterCharge_exp_sq_p15_std_10v_st1000_nt200, "exp, cluster", "pe");
//     legend->AddEntry(h_clusterCharge_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, cluster", "pe");
//     legend->Draw();

//     // TLatex title;
//     // title.SetTextAlign(12);
//     // title.SetTextSize(0.05);
//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 STD 10V");

//     canvas->SaveAs("./plot/clusterCharge_sim_exp_std.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(20);
//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(20);
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);

//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed+1);
//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed+1);
//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(20);
//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed+1);
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed+1);
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(20);
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);

//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);

//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 2000);
//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 2000);

//     h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterNeighborsChargeSum_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, neighbor sum", "pe");
//     legend->AddEntry(h_clusterNeighborsChargeSum_exp_sq_p15_std_10v_st1000_nt200, "exp, neighbor sum", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 STD 10V");
//     canvas->SaveAs("./plot/clusterNeighborsChargeSum_sim_exp_std.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterNeighborsChargeSum_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, neighbor sum", "pe");
//     legend->AddEntry(h_clusterNeighborsChargeSum_exp_sq_p15_gap_10v_st1000_nt200, "exp, neighbor sum", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.6, 0.85, "SQ P15 GAP 10V");
//     canvas->SaveAs("./plot/clusterNeighborsChargeSum_sim_exp_gap.pdf");


//     // Seed Charge for each cluster size
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size1_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size1_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 GAP 10V cs1");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs1_sim_exp_gap.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size2_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size2_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 GAP 10V cs2");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs2_sim_exp_gap.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size3_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size3_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 GAP 10V cs3");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs3_sim_exp_gap.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size4_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size4_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 GAP 10V cs4");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs4_sim_exp_gap.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 GAP 10V cs5");

//     legend->AddEntry(h_clusterSeedCharge_size5_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size5_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs5_sim_exp_gap.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size6_exp_sq_p15_gap_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size6_sim_sq_p15_gap_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 GAP 10V cs6");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs6_sim_exp_gap.pdf");

//  // ===================================================================================
//  // ===================================================================================
//     // Seed Charge for each cluster size
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size1_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size1_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 STD 10V cs1");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs1_sim_exp_std.pdf");

//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size2_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size2_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 STD 10V cs2");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs2_sim_exp_std.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size3_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size3_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 STD 10V cs3");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs3_sim_exp_std.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 STD 10V cs4");

//     legend->AddEntry(h_clusterSeedCharge_size4_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size4_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs4_sim_exp_std.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());

//     h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size5_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size5_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 STD 10V cs5");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs5_sim_exp_std.pdf");
//     canvas->Clear();
//     legend->Clear();

//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerColor(kRed);
//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetLineColor(kRed);
//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->SetMarkerStyle(24);

//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->SetMarkerColor(kBlue);
//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->SetLineColor(kBlue);
//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->SetMarkerSize(0.8);
//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->SetMarkerStyle(24);

//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->Rebin(10);
//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetXaxis()->SetRangeUser(0, 4000);

//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Scale(1.0/h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->GetMaximum());
//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->Scale(1.0/h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->GetMaximum());

//     // legend->AddEntry(h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     // legend->AddEntry(h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     // legend->Draw();

//     h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um->Draw("PE");
//     h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200->Draw("samePE");

//     legend->AddEntry(h_clusterSeedCharge_size6_exp_sq_p15_std_10v_st1000_nt200, "exp, seed", "pe");
//     legend->AddEntry(h_clusterSeedCharge_size6_sim_sq_p15_std_10v_st0_nt0_n0e_1um, "sim, seed", "pe");
//     legend->Draw();
    

//     title.DrawLatexNDC(0.5, 0.85, "SQ P15 STD 10V cs6");

//     canvas->SaveAs("./plot/seed_charge_cs/clusterSeedCharge_cs6_sim_exp_std.pdf");

// } 

// void plot_Check::run_PlotCheck() {
    
// }
