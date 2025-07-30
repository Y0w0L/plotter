#include "plot_ExperimentData.h"

plot_ExperimentData::plot_ExperimentData() {
    LOG_STATUS.source("plot_ExperimentData::plot_ExperimentData") << "Plot_ExperimentData object is created.";
    PIXEL_PITCH_ = {"22p5"};
    CHIP_TYPE_ = {"std", "gap"};
    //CHIP_TYPE_ = {"blk"};
    VOLTAGE_ = {"10", "7", "4"};
    //VOLTAGE_ = {"10"};
    //SEED_THRESHOLD_ = {"400"};
    //SEED_THRESHOLD_ = {"300"};
    SEED_THRESHOLD_ = {"500"};
    //NEIGHBOR_THRESHOLD_ = {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"}; 
    //NEIGHBOR_THRESHOLD_ = {"50","60" ,"70", "80", "90", "100"}; 
    NEIGHBOR_THRESHOLD_ = {"50", "80", "100", "200", "300", "350", "400", "450", "500"};
    //NEIGHBOR_THRESHOLD_ = {"60"};
    TIME_ = plot_histogram::currentDateTime();
}

std::vector<int> plot_ExperimentData::GetMyColors(int n_colors) {
    gStyle->SetPalette(kRainBow);

    std::vector<int> colors;
    for(int i=0; i<n_colors; ++i) {
        int color_index_in_palette = int((i + 0.5) * (gStyle->GetNumberOfColors() / double(n_colors)));
        colors.push_back(TColor::GetColorPalette(color_index_in_palette));
    }
    return colors;
}

void plot_ExperimentData::set_2DSURFStyle(TCanvas* canvas, TH2D* hist) {
    canvas->Clear();
    canvas->SetTopMargin(0.062);
    canvas->SetBottomMargin(0.1);
    canvas->SetLeftMargin(0.14);
    canvas->SetRightMargin(0.04);

    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetZaxis()->SetTitleOffset(0.8);

    hist->Draw("SURF1");
}

void plot_ExperimentData::set_GraphStyle(TGraph* graph) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetLineWidth(1);
}

void plot_ExperimentData::set_1DStyle(TH1D* hist) {
    hist->SetStats(0);
    hist->GetXaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->SetTitleFont(42, "XYZ");

    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(20);
}

TH2D* plot_ExperimentData::convert_toTH2D(TProfile2D* profile2D) {
    if(!profile2D) {
        return nullptr;
    }

    // get axis information from TProfile2D
    int nbinsX  = profile2D->GetNbinsX();
    double xMin = profile2D->GetXaxis()->GetXmin();
    double xMax = profile2D->GetXaxis()->GetXmax();
    // double xMin = -nbinsX / 2;
    // double xMax = nbinsX / 2;

    int nbinsY  = profile2D->GetNbinsY();
    double yMin = profile2D->GetYaxis()->GetXmin();
    double yMax = profile2D->GetYaxis()->GetXmax();
    // double yMin = -nbinsY / 2;
    // double yMax = nbinsY / 2;

    // make new TH2D
    TH2D* h2d = new TH2D(profile2D->GetName(), profile2D->GetTitle(), nbinsX, xMin, xMax, nbinsY, yMin, yMax);

    double content;
    double error;
    // copy contents and fill
    for(int i=1; i<=nbinsX; i++) {
        for(int j=1; j<=nbinsY; j++) {
            content = profile2D->GetBinContent(i,j);
            error = profile2D->GetBinError(i,j);

            h2d->SetBinContent(i,j,content);
            h2d->SetBinError(i,j,error);
        }
    }

    h2d->SetStats(0);
    h2d->GetXaxis()->SetTitleOffset(0.8);
    h2d->GetYaxis()->SetTitleOffset(0.7);
    h2d->GetZaxis()->SetTitleOffset(0.6);
    h2d->GetXaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetZaxis()->SetTitleSize(0.06);
    h2d->GetXaxis()->SetLabelSize(0.04);
    h2d->GetYaxis()->SetLabelSize(0.04);
    h2d->GetYaxis()->SetLabelSize(0.04);

    h2d->SetContour(20);

    return h2d;
}

TH1D* plot_ExperimentData::get_TH1D(std::string filename) {
    // TFile* file = TFile::Open(filename);
    // TH1D* hist = (TH1D*)file->Get("AnalysisCE65/CE65_3/local_residuals/residualsX");
    // return hist;
}

void plot_ExperimentData::run_NoiseScan() {
    LOG_STATUS.source("plot_ExperimentData::run_NoiseScan") << "Start run for Noise scan data.";

    std::string output_file_name = "/home/towa/alice3/plotter/plot/experimentData_NoiseScan.root";
    std::string data_dir_path = "/home/towa/alice3/hist/kek202412/";
    TFile* output = TFile::Open(output_file_name.c_str(), "RECREATE");

    TDirectory* pixelCharge = output->mkdir("pixelCharge");

    std::vector<TDirectory*> outputDir = {pixelCharge};

    TFile* inputROOTFile;
    TH1D* h_pixelCharge = nullptr;
    std::vector<TH1D*> v_pixelCharge = {};

    TCanvas* canvas = new TCanvas("canvas","canvas",800,600);
    TLatex title;
    TLatex condition;
    title.SetTextSize(0.04);
    title.SetTextFont(62);
    condition.SetTextSize(0.03);
    condition.SetTextFont(62);

    std::string chip_variation;
    std::string chip_variation_text;

    for(int i=0; i<PIXEL_PITCH_.size(); i++) {
        for(int j=0; j<CHIP_TYPE_.size(); j++) {
            canvas->Clear();
            canvas->SetTopMargin(0.062);
            canvas->SetBottomMargin(0.14);
            canvas->SetLeftMargin(0.13);
            canvas->SetRightMargin(0.07);
            for(int k=0; k<VOLTAGE_.size(); k++) {
                inputROOTFile = TFile::Open(Form("%skek202412_%s_%s_%sV_SeedThd0e_NeighborThd0e_noise.root", data_dir_path.c_str(), PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str()));
                chip_variation = Form("%s_%s_%sV_SeedThd0e_NeighborThd0e", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str());
                chip_variation_text = Form("p%s/%s/%sV/SeedThd0e/NeighborThd0e", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str());
                gStyle->SetPalette(kViridis);

                h_pixelCharge = (TH1D*)inputROOTFile->Get("EventLoaderEUDAQ2/CE65_3/hPixelRawValues");
                h_pixelCharge->SetTitle(";charge [ADC];counts");
                plot_ExperimentData::set_1DStyle(h_pixelCharge);
                if(k==0) {
                    h_pixelCharge->GetXaxis()->SetRangeUser(0,120);
                    //h_pixelCharge->GetYaxis()->SetRangeUser(0,2100000);
                    h_pixelCharge->SetMarkerColor(kBlack);
                    h_pixelCharge->SetLineColor(kBlack);
                    h_pixelCharge->Draw("PE");
                    v_pixelCharge.push_back(h_pixelCharge);
                } else {
                    if(k==1) {
                        h_pixelCharge->SetMarkerColor(kRed);
                        h_pixelCharge->SetLineColor(kRed);
                    }
                    if(k==2) {
                        h_pixelCharge->SetMarkerColor(kBlue);
                        h_pixelCharge->SetLineColor(kBlue);
                    }
                    h_pixelCharge->Draw("samePE");
                    v_pixelCharge.push_back(h_pixelCharge);
                }
                //inputROOTFile->Close();
            } // VOLTAGE_

            title.DrawLatexNDC(0.70, 0.89, "Per-Pixel Charge");
            condition.DrawLatexNDC(0.70, 0.85, "w/o beam");
            condition.DrawLatexNDC(0.70, 0.79, Form("%s/%s", CHIP_TYPE_[j].c_str(), PIXEL_PITCH_[i].c_str()));
            condition.DrawLatexNDC(0.70, 0.82, Form("Plotted on %s", TIME_.c_str()));

            TLegend* legend = new TLegend(0.75, 0.60, 0.85, 0.78);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.04);
            legend->SetBorderSize(0);
            for(int k=0; k<VOLTAGE_.size(); k++) {
                legend->AddEntry(v_pixelCharge[k], Form("%sV", VOLTAGE_[k].c_str()), "PE");
            }
            legend->Draw();

            pixelCharge->cd();
            canvas->Write(Form("perPixelCharge_ce65_%s", chip_variation.c_str()));
            canvas->SaveAs(Form("plot/perPixelCharge_ce65_%s.pdf", chip_variation.c_str()));
        } // CHIP_TYPE_
    } // PIXEL_PITCH_

    output->Close();
}

void plot_ExperimentData::run_Analysis() {
    LOG_STATUS.source("plot_ExperimentData::run_Analysis") << "Start run for beamtest analysis.";

    std::string output_file_name = "/home/towa/alice3/plotter/plot/experimentData_analysis.root";
    std::string data_dir_path = "/home/towa/alice3/hist/kek202412/";
    TFile* output = TFile::Open(output_file_name.c_str(), "RECREATE");

    TDirectory* clusterCharge = output->mkdir("cluster_charge");
    TDirectory* seedCharge = output->mkdir("seed_charge");
    TDirectory* clusterSize = output->mkdir("cluster_size");
    TDirectory* residual = output->mkdir("residual");
    TDirectory* efficiency = output->mkdir("efficiency");
    TDirectory* tgraph = output->mkdir("tgraph");

    std::vector<int> my_colors = plot_ExperimentData::GetMyColors(NEIGHBOR_THRESHOLD_.size());

    std::vector<TDirectory*> outputDir = {clusterCharge, seedCharge, clusterSize, residual, efficiency};

    TFile* inputROOTFile;
    TH1D* hClusterCharge = nullptr;
    TH1D* hSeedCharge = nullptr;
    TH1D* hClusterSize = nullptr;
    TH1D* hResidualX = nullptr;
    TH1D* hResidualY = nullptr;
    //TH1D* hResidualR = nullptr;
    //TH1D* hEfficiency = nullptr;
    TF1* fResidualX = nullptr;

    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    TLatex title;
    TLatex condition;
    title.SetTextSize(0.04);
    title.SetTextFont(62);
    condition.SetTextSize(0.03);
    condition.SetTextFont(62);

    gStyle->SetPalette(55);

    std::string chip_variation;
    std::string chip_variation_text;

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vMeanClusterSize;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vResolutionX;
    vResolutionX.resize(PIXEL_PITCH_.size());
    vMeanClusterSize.resize(PIXEL_PITCH_.size());

    for(size_t i=0; i < PIXEL_PITCH_.size(); i++) {
        vResolutionX[i].resize(CHIP_TYPE_.size());
        vMeanClusterSize[i].resize(CHIP_TYPE_.size());

        for(size_t j=0; j < CHIP_TYPE_.size(); j++) {
            vResolutionX[i][j].resize(VOLTAGE_.size());
            vMeanClusterSize[i][j].resize(VOLTAGE_.size());

            for(size_t k=0; k < VOLTAGE_.size(); k++) {
                vResolutionX[i][j][k].resize(SEED_THRESHOLD_.size());
                vMeanClusterSize[i][j][k].resize(SEED_THRESHOLD_.size());
            }
        }
    }

    for(int i=0; i<PIXEL_PITCH_.size(); i++) {
        for(int j=0; j<CHIP_TYPE_.size(); j++) {
            for(int k=0; k<VOLTAGE_.size(); k++) {
                for(int l=0; l<SEED_THRESHOLD_.size(); l++) {
                    canvas->Clear();
                    canvas->SetTopMargin(0.062);
                    canvas->SetBottomMargin(0.14);
                    canvas->SetLeftMargin(0.13);
                    canvas->SetRightMargin(0.07);

                    TLegend* legend = new TLegend(0.65, 0.50, 0.85, 0.78);
                    legend->SetFillStyle(0);
                    legend->SetTextSize(0.04);
                    legend->SetBorderSize(0);

                    chip_variation = Form("%s_%s_%sV_SeedThd%se", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str());
                    chip_variation_text = Form("p%s/%s/%sV/SeedThd%s", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str());
                        

                    for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                        inputROOTFile = TFile::Open(Form("%skek202412_%s_%s_%sV_SeedThd%se_NeighborThd%se.root", data_dir_path.c_str(), PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str()));
                        //chip_variation = Form("%s_%s_%sV_SeedThd%se_NeighborThd%se", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str());
                        //chip_variation_text = Form("p%s/%s/%sV/SeedThd%s/NeighborThd%s", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str());
                        //gStyle->SetPalette(kViridis);

                        hClusterCharge = (TH1D*)inputROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
                        hSeedCharge = (TH1D*)inputROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterSeedCharge");
                        hClusterCharge->SetTitle(";charge [ADC];counts");
                        hClusterCharge->Rebin(10);
                        plot_ExperimentData::set_1DStyle(hClusterCharge);
                        hClusterCharge->SetMarkerColorAlpha(my_colors[n], 1);
                        hClusterCharge->SetLineColorAlpha(my_colors[n], 1);
                        legend->AddEntry(hClusterCharge, Form("Thd%s", NEIGHBOR_THRESHOLD_[n].c_str()), "PE");
                        if(n==0) {
                            hClusterCharge->GetXaxis()->SetRangeUser(50, 5000);
                            //hClusterCharge->GetYaxis()->SetRangeUser(0, 1400);
                            // hClusterCharge->SetMarkerColor(kRed);
                            // hClusterCharge->SetLineColor(kRed);
                            hClusterCharge->Draw("PE");
                            hClusterCharge->Draw("same C");
                        } else {
                            hClusterCharge->Draw("samePE");
                            hClusterCharge->Draw("same C");
                        }
                    } // NEIGHBOR_THRESHOLD_

                    title.DrawLatexNDC(0.60, 0.89, "Cluster Charge");
                    condition.DrawLatexNDC(0.60, 0.85, "Electron 3GeV/c @KEK");
                    condition.DrawLatexNDC(0.60, 0.82, chip_variation_text.c_str());
                    condition.DrawLatexNDC(0.60, 0.79, Form("Plotted on %s", TIME_.c_str()));
                    legend->Draw();

                    clusterCharge->cd();
                    canvas->Write(Form("clusterCharge_ce65_%s", chip_variation.c_str()));
                    canvas->SaveAs(Form("plot/clusterCharge_ce65_%s.pdf", chip_variation.c_str()));

                    // Position Resolution with Cluster Size
                    legend->Clear();
                    canvas->Clear();
                    canvas->SetTopMargin(0.062);
                    canvas->SetBottomMargin(0.14);
                    canvas->SetLeftMargin(0.13);
                    canvas->SetRightMargin(0.07);

                    //TLegend* legend = new TLegend(0.65, 0.50, 0.85, 0.78);
                    legend->SetFillStyle(0);
                    legend->SetTextSize(0.04);
                    legend->SetBorderSize(0);

                    for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                        inputROOTFile = TFile::Open(Form("%skek202412_%s_%s_%sV_SeedThd%se_NeighborThd%se.root", data_dir_path.c_str(), PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str()));
                        
                        hResidualX = (TH1D*)inputROOTFile->Get("AnalysisCE65/CE65_3/local_residuals/residualsX");
                        hClusterSize = (TH1D*)inputROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterSize");
                        hResidualX->SetTitle(";x_{track} - x_{hit} [um];counts");
                        hResidualX->Rebin(10);
                        hResidualX->Scale(1/hResidualX->GetMaximum());
                        plot_ExperimentData::set_1DStyle(hResidualX);
                        hResidualX->SetMarkerColorAlpha(my_colors[n], 1);
                        hResidualX->SetLineColorAlpha(my_colors[n], 1);
                        legend->AddEntry(hResidualX, Form("Thd%s", NEIGHBOR_THRESHOLD_[n].c_str()), "PE");

                        if(n==0) {
                            hResidualX->GetXaxis()->SetRangeUser(-50, 50);
                            hResidualX->GetYaxis()->SetRangeUser(0, 1.1);
                            // hResidualX->SetMarkerColor(kRed);
                            // hResidualX->SetLineColor(kRed);
                            hResidualX->Draw("PE");
                        } else {
                            hResidualX->Draw("same PE");
                        }

                        fResidualX = new TF1("fResidualX", "gaus", -50, 50);
                        fResidualX->SetLineColor(my_colors[n]);
                        //fResidualX->SetLineColor(kRed);
                        hResidualX->Fit(fResidualX, "RQ");

                        vResolutionX[i][j][k][l].push_back(fResidualX->GetParameter(2));
                        vMeanClusterSize[i][j][k][l].push_back(hClusterSize->GetMean());
                    }// NEIGHBOR_THRESHOLD_

                    title.DrawLatexNDC(0.60, 0.89, "Residual (X_{track} - X_{hit})");
                    condition.DrawLatexNDC(0.60, 0.85, "Electron 3GeV/c @KEK");
                    condition.DrawLatexNDC(0.60, 0.82, chip_variation_text.c_str());
                    condition.DrawLatexNDC(0.60, 0.79, Form("Plotted on %s", TIME_.c_str()));
                    legend->Draw();

                    residual->cd();
                    canvas->Write(Form("residualX_ce65_%s", chip_variation.c_str()));
                    canvas->SaveAs(Form("plot/residualX_ce65_%s.pdf", chip_variation.c_str()));
                } // SEED_THRESHOLD_
            } // VOLTAGE_
        } // CHIP_TYPE_
    } // PIXEL_PTICH_

    bool tgraph_each_chip = false;
    if(tgraph_each_chip) {
        for(int i=0; i<PIXEL_PITCH_.size(); i++) {
            for(int j=0; j<CHIP_TYPE_.size(); j++) {
                for(int l=0; l<SEED_THRESHOLD_.size(); l++) {

                    // --- 1. 分解能 (Resolution) のプロット ---
                    canvas->Clear();
                    canvas->SetGrid();

                    TMultiGraph* mg_resolution = new TMultiGraph();
                    TLegend* legend_res = new TLegend(0.65, 0.70, 0.88, 0.88);
                    legend_res->SetFillStyle(0);
                    legend_res->SetBorderSize(0);

                    // 電圧(k)ごとにグラフを作成し、MultiGraphに追加
                    for(int k=0; k<VOLTAGE_.size(); k++) {
                        TGraph* gr_resolution = new TGraph();

                        // データをセット
                        for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                            double x_val = std::stod(NEIGHBOR_THRESHOLD_[n]);
                            double y_resolution = vResolutionX[i][j][k][l][n];
                            gr_resolution->SetPoint(n, x_val, y_resolution);
                        }

                        // 電圧ごとに色とスタイルを設定
                        gr_resolution->SetMarkerStyle(20 + k);
                        gr_resolution->SetMarkerColor(kBlue - 3*k);
                        gr_resolution->SetLineColor(kBlue - 3*k);
                        gr_resolution->SetLineWidth(2);

                        mg_resolution->Add(gr_resolution, "PL");
                        legend_res->AddEntry(gr_resolution, Form("%sV", VOLTAGE_[k].c_str()), "pl");
                    }

                    // MultiGraphを描画
                    mg_resolution->Draw("A");
                    mg_resolution->SetTitle(";Neighbor Threshold [e];Position Resolution #sigma [#mum]");

                    // 軸のスタイルを設定
                    mg_resolution->GetXaxis()->SetTitleSize(0.05);
                    mg_resolution->GetXaxis()->SetLabelSize(0.04);
                    mg_resolution->GetYaxis()->SetTitleSize(0.05);
                    mg_resolution->GetYaxis()->SetLabelSize(0.04);

                    legend_res->Draw();
                    chip_variation_text = Form("p%s / %s / SeedThd %se", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), SEED_THRESHOLD_[l].c_str());
                    title.DrawLatexNDC(0.15, 0.91, chip_variation_text.c_str());

                    canvas->SaveAs(Form("plot/Resolution_vs_NeighborThd_CompV_%s.pdf", Form("%s_%s_Seed%s",PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), SEED_THRESHOLD_[l].c_str())));
                    tgraph->cd();
                    canvas->Write(Form("Resolution_vs_NeighborThd_CompV_%s", Form("%s_%s_Seed%s",PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), SEED_THRESHOLD_[l].c_str())));


                    // --- 2. クラスターサイズ (Cluster Size) のプロット ---
                    canvas->Clear();
                    canvas->SetGrid();

                    TMultiGraph* mg_clustersize = new TMultiGraph();
                    TLegend* legend_cs = new TLegend(0.65, 0.70, 0.88, 0.88);
                    legend_cs->SetFillStyle(0);
                    legend_cs->SetBorderSize(0);

                    for(int k=0; k<VOLTAGE_.size(); k++) {
                        TGraph* gr_clustersize = new TGraph();

                        for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                            double x_val = std::stod(NEIGHBOR_THRESHOLD_[n]);
                            double y_clustersize = vMeanClusterSize[i][j][k][l][n];
                            gr_clustersize->SetPoint(n, x_val, y_clustersize);
                        }

                        gr_clustersize->SetMarkerStyle(20 + k);
                        gr_clustersize->SetMarkerColor(kRed - 3*k);
                        gr_clustersize->SetLineColor(kRed - 3*k);
                        gr_clustersize->SetLineWidth(2);

                        mg_clustersize->Add(gr_clustersize, "PL");
                        legend_cs->AddEntry(gr_clustersize, Form("%sV", VOLTAGE_[k].c_str()), "pl");
                    }

                    mg_clustersize->Draw("A");
                    mg_clustersize->SetTitle(";Neighbor Threshold [e];Mean Cluster Size");

                    mg_clustersize->GetXaxis()->SetTitleSize(0.05);
                    mg_clustersize->GetXaxis()->SetLabelSize(0.04);
                    mg_clustersize->GetYaxis()->SetTitleSize(0.05);
                    mg_clustersize->GetYaxis()->SetLabelSize(0.04);

                    legend_cs->Draw();
                    title.DrawLatexNDC(0.15, 0.91, chip_variation_text.c_str());

                    canvas->SaveAs(Form("plot/ClusterSize_vs_NeighborThd_CompV_%s.pdf", Form("%s_%s_Seed%s",PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), SEED_THRESHOLD_[l].c_str())));
                    tgraph->cd();
                    canvas->Write(Form("ClusterSize_vs_NeighborThd_CompV_%s", Form("%s_%s_Seed%s",PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), SEED_THRESHOLD_[l].c_str())));


                }
            }
        }
    }

    // bool residual_plots = true;
    // if(residual_plots) {
    //     canvas->Clear();
    //     canvas->SetTopMargin(0.062);
    //     canvas->SetBottomMargin(0.14);
    //     canvas->SetLeftMargin(0.13);
    //     canvas->SetRightMargin(0.07);

    //     TH1D* h_std_4v_60 = get_TH1D("")

    // } // residual_plots

    bool tgraph_all = false;
    if(tgraph_all) {
                // --- 1. 分解能 (Resolution) のプロット ---
        for(int i=0; i<PIXEL_PITCH_.size(); i++) {
            for(int l=0; l<SEED_THRESHOLD_.size(); l++) {
                canvas->Clear();
                canvas->SetTopMargin(0.062);
                canvas->SetBottomMargin(0.14);
                canvas->SetLeftMargin(0.13);
                canvas->SetRightMargin(0.07);
                canvas->SetGrid();
                    
                TMultiGraph* mg_resolution = new TMultiGraph();
                TLegend* legend_res = new TLegend(0.55, 0.65, 0.88, 0.88); // 凡例のサイズを調整
                legend_res->SetFillStyle(0);
                legend_res->SetBorderSize(0);
                    
                // j (チップ) と k (電圧) の両方でループ
                for(int j=0; j<CHIP_TYPE_.size(); j++) {
                    for(int k=0; k<VOLTAGE_.size(); k++) {
                        TGraph* gr_resolution = new TGraph();
                    
                        for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                            double x_val = std::stod(NEIGHBOR_THRESHOLD_[n]);
                            double y_resolution = vResolutionX[i][j][k][l][n];
                            gr_resolution->SetPoint(n, x_val, y_resolution);
                        }
                    
                        // 色で電圧、線の種類でチップを表現
                        gr_resolution->SetMarkerStyle(20 + k + j*VOLTAGE_.size());
                        gr_resolution->SetMarkerColor(kBlue - 3*k);
                        gr_resolution->SetLineColor(kBlue - 3*k);
                        gr_resolution->SetLineStyle(j + 1); // 1: 実線, 2: 破線
                        gr_resolution->SetLineWidth(1.5);
                    
                        mg_resolution->Add(gr_resolution, "PL");
                        legend_res->AddEntry(gr_resolution, Form("%s, %sV", CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str()), "pl");
                    }
                }
            
                mg_resolution->Draw("A");
                mg_resolution->SetTitle(";threshold [ADC];resolution in x [um]");
                mg_resolution->GetYaxis()->SetRangeUser(0, 10);
                mg_resolution->GetXaxis()->SetLimits(40,510);
                //mg_resolution->GetYaxis()->SetLimits(0,10);
                mg_resolution->GetXaxis()->SetTitleSize(0.06);
                mg_resolution->GetXaxis()->SetLabelSize(0.04);
                mg_resolution->GetXaxis()->SetTitleOffset(0.6);
                mg_resolution->GetYaxis()->SetTitleSize(0.06);
                mg_resolution->GetYaxis()->SetLabelSize(0.04);
                mg_resolution->GetYaxis()->SetTitleOffset(0.7);
            
                legend_res->Draw();
                chip_variation_text = Form("p%s/SeedThd%sADC", PIXEL_PITCH_[i].c_str(), SEED_THRESHOLD_[l].c_str());
                title.DrawLatexNDC(0.15, 0.91, chip_variation_text.c_str());
            
                canvas->SaveAs(Form("plot/Resolution_CompAll_%s.pdf", Form("%s_Seed%s",PIXEL_PITCH_[i].c_str(), SEED_THRESHOLD_[l].c_str())));
            
                        // --- 2. クラスターサイズ (Cluster Size) のプロット ---
                canvas->Clear();
                canvas->SetTopMargin(0.062);
                canvas->SetBottomMargin(0.14);
                canvas->SetLeftMargin(0.13);
                canvas->SetRightMargin(0.07);
                canvas->SetGrid();
            
                TMultiGraph* mg_clustersize = new TMultiGraph();
                TLegend* legend_cs = new TLegend(0.55, 0.65, 0.88, 0.88);
                legend_cs->SetFillStyle(0);
                legend_cs->SetBorderSize(0);
            
                for(int j=0; j<CHIP_TYPE_.size(); j++) {
                    for(int k=0; k<VOLTAGE_.size(); k++) {
                        TGraph* gr_clustersize = new TGraph();
                    
                        for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                            double x_val = std::stod(NEIGHBOR_THRESHOLD_[n]);
                            double y_clustersize = vMeanClusterSize[i][j][k][l][n];
                            gr_clustersize->SetPoint(n, x_val, y_clustersize);
                        }
                    
                        gr_clustersize->SetMarkerStyle(20 + k + j*VOLTAGE_.size());
                        gr_clustersize->SetMarkerColor(kRed - 3*k);
                        gr_clustersize->SetLineColor(kRed - 3*k);
                        gr_clustersize->SetLineStyle(j + 1); // 1: 実線, 2: 破線
                        gr_clustersize->SetLineWidth(1.5);
                    
                        mg_clustersize->Add(gr_clustersize, "PL");
                        legend_cs->AddEntry(gr_clustersize, Form("%s, %sV", CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str()), "pl");
                    }
                }
            
                mg_clustersize->Draw("A");
                mg_clustersize->SetTitle(";threshold [ADC];mean cluster size");
                mg_clustersize->GetXaxis()->SetLimits(40,510);
                mg_clustersize->GetYaxis()->SetRangeUser(0, 6);
                mg_clustersize->GetXaxis()->SetTitleSize(0.06);
                mg_clustersize->GetXaxis()->SetLabelSize(0.04);
                mg_clustersize->GetXaxis()->SetTitleOffset(0.6);
                mg_clustersize->GetYaxis()->SetTitleSize(0.06);
                mg_clustersize->GetYaxis()->SetLabelSize(0.04);
                mg_clustersize->GetYaxis()->SetTitleOffset(0.7);
            
                legend_cs->Draw();
                title.DrawLatexNDC(0.15, 0.91, chip_variation_text.c_str()); // chip_variation_textは前のプロットから流用
            
                canvas->SaveAs(Form("plot/ClusterSize_CompAll_%s.pdf", Form("%s_Seed%s",PIXEL_PITCH_[i].c_str(), SEED_THRESHOLD_[l].c_str())));
            
                    // --- 3. 2軸比較プロット (左: Resolution, 右: Cluster Size) ---
                canvas->Clear();
                canvas->SetTopMargin(0.052);
                canvas->SetBottomMargin(0.145);
                canvas->SetLeftMargin(0.13);
                canvas->SetRightMargin(0.13);
            
                TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
                // pad1->SetBottomMargin(0.14);
                // pad1->SetRightMargin(0.16);
                // pad1->SetLeftMargin(0.16);
                pad1->SetGrid();
                pad1->Draw();
            
                canvas->cd();
                TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
                pad2->SetFillStyle(0);
                pad2->SetFrameFillStyle(0);
                pad2->Draw();
            
                TMultiGraph* mg_res_comp = new TMultiGraph();
                TLegend* legend_comp = new TLegend(0.70, 0.4, 0.93, 0.63);
                legend_comp->SetFillStyle(0);
                legend_comp->SetBorderSize(0);
                TLegend* legend_comp_clsize = new TLegend(0.65, 0.4, 0.88, 0.63);
                legend_comp_clsize->SetFillStyle(0);
                legend_comp_clsize->SetBorderSize(0);
            
                std::vector<TGraph*> v_gr_cs;
            
                for(int j=0; j<CHIP_TYPE_.size(); j++) {
                    for(int k=0; k<VOLTAGE_.size(); k++) {
                        TGraph* gr_res = new TGraph();
                        TGraph* gr_cs = new TGraph();
                    
                        for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                            double x_val = std::stod(NEIGHBOR_THRESHOLD_[n]);
                            double y_res = vResolutionX[i][j][k][l][n];
                            double y_cs = vMeanClusterSize[i][j][k][l][n];
                            gr_res->SetPoint(n, x_val, y_res);
                            gr_cs->SetPoint(n, x_val, y_cs);
                        }
                    
                        // Resolution (左軸) のスタイル
                        gr_res->SetMarkerStyle(20);
                        gr_res->SetMarkerColor(kBlue - 3*k);
                        gr_res->SetLineColor(kBlue - 3*k);
                        gr_res->SetLineStyle(j + 1);
                        gr_res->SetLineWidth(2);
                    
                        // Cluster Size (右軸) のスタイル
                        gr_cs->SetMarkerStyle(21);
                        gr_cs->SetMarkerColor(kRed - 3*k);
                        gr_cs->SetLineColor(kRed - 3*k);
                        gr_cs->SetLineStyle(j + 1);
                        gr_cs->SetLineWidth(2);
                    
                        mg_res_comp->Add(gr_res, "PL");
                        v_gr_cs.push_back(gr_cs);
                        legend_comp->AddEntry(gr_cs, Form("%s, %sV", CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str()), "pl");
                        legend_comp_clsize->AddEntry(gr_res, "", "pl");
                    }
                }
            
                // --- 左軸 (Resolution) の描画 ---
                pad1->cd();
                mg_res_comp->Draw("A");
                mg_res_comp->SetTitle(";threshold [ADC];resolution in x [um]");
                mg_res_comp->GetXaxis()->SetLimits(40,510);
                mg_res_comp->GetYaxis()->SetLimits(3.5,9.5);
                mg_res_comp->GetXaxis()->SetTitleSize(0.06);
                mg_res_comp->GetXaxis()->SetLabelSize(0.04);
                mg_res_comp->GetXaxis()->SetTitleOffset(0.7);
                mg_res_comp->GetYaxis()->SetTitleSize(0.06);
                mg_res_comp->GetYaxis()->SetLabelSize(0.04);
                mg_res_comp->GetYaxis()->SetTitleOffset(0.6);
                //mg_res_comp->GetYaxis()->SetTitleColor(kBlue);
                //mg_res_comp->GetYaxis()->SetLabelColor(kBlue);
            
                // --- 右軸 (Cluster Size) の描画 ---
                pad2->cd();
                double cs_min = 1e9, cs_max = -1e9;
                for(auto gr : v_gr_cs) {
                    if (gr->GetYaxis()->GetXmin() < cs_min) cs_min = gr->GetYaxis()->GetXmin();
                    if (gr->GetYaxis()->GetXmax() > cs_max) cs_max = gr->GetYaxis()->GetXmax();
                }
            
                for(size_t ig=0; ig<v_gr_cs.size(); ++ig) {
                    if (ig==0) {
                        v_gr_cs[ig]->Draw("APL"); 
                        v_gr_cs[ig]->GetXaxis()->SetRangeUser(40, 510);
                        v_gr_cs[ig]->GetYaxis()->SetRangeUser(0.5, 6.5);
                        v_gr_cs[ig]->GetXaxis()->SetLabelSize(0);
                        v_gr_cs[ig]->GetXaxis()->SetTitleSize(0);
                        v_gr_cs[ig]->GetYaxis()->SetLabelSize(0);
                        v_gr_cs[ig]->GetYaxis()->SetTitleSize(0);
                    } else {
                        v_gr_cs[ig]->Draw("PL same");
                    }
                }
            
                pad2->Update();
                TGaxis *axis_cs_comp = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0.5, 6.5, 510, "+L");
                axis_cs_comp->SetTitle("mean cluster size");
                // axis_cs_comp->SetTitleColor(kRed);
                // axis_cs_comp->SetLineColor(kRed);
                // axis_cs_comp->SetLabelColor(kRed);
                axis_cs_comp->SetTitleSize(0.06);
                axis_cs_comp->SetLabelSize(0.04);
                axis_cs_comp->SetTitleOffset(0.6);
                axis_cs_comp->SetTitleFont(42);
                axis_cs_comp->SetLabelFont(42);
                axis_cs_comp->Draw();
            
                canvas->cd();
                legend_comp->Draw();
                legend_comp_clsize->Draw();
                //title.DrawLatexNDC(0.15, 0.91, chip_variation_text.c_str()); // chip_variation_textは流用
                condition.DrawLatexNDC(0.12, 0.87, "e^{-} 3GeV/c @KEK(Dec. 2024)");
                condition.DrawLatexNDC(0.12, 0.84, chip_variation_text.c_str());
                condition.DrawLatexNDC(0.65,0.87, Form("Plotted on %s", TIME_.c_str()));
            
                canvas->SaveAs(Form("plot/DualAxis_CompAll_%s.pdf", Form("%s_Seed%s",PIXEL_PITCH_[i].c_str(), SEED_THRESHOLD_[l].c_str())));
        
            }
        }
    }


    // writing about blk cluster charge
    bool blk_check = false;
    if(blk_check) {
        canvas->Clear();
        canvas->SetTopMargin(0.062);
        canvas->SetBottomMargin(0.14);
        canvas->SetLeftMargin(0.13);
        canvas->SetRightMargin(0.07);

        std::vector<TH1D*> h_blks;

        TFile* file_std = TFile::Open(Form("%skek202412_22p5_std_10V_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str()));
        TFile* file_gap = TFile::Open(Form("%skek202412_22p5_blk_10V_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str()));
        TH1D* h_std = (TH1D*)file_std->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
        TH1D* h_gap = (TH1D*)file_gap->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");

        h_std->SetTitle(";charge [ADC];#counts");
        plot_ExperimentData::set_1DStyle(h_std);
        plot_ExperimentData::set_1DStyle(h_gap);
        h_std->SetMarkerColor(kGray+1);
        h_std->SetLineColor(kGray+1);
        h_gap->SetMarkerColor(kGray+2);
        h_gap->SetLineColor(kGray+2);

        h_std->Rebin(10);
        h_gap->Rebin(10);
        h_std->Scale(1/h_std->GetMaximum());
        h_gap->Scale(1/h_gap->GetMaximum());

        h_std->GetXaxis()->SetRangeUser(50, 4000);
        h_std->GetYaxis()->SetRangeUser(0, 1.1);

        h_std->Draw("PE");
        h_std->Draw("same C");

        h_gap->Draw("same PE");
        h_gap->Draw("same C");

        TFile* file_blk_4v = TFile::Open(Form("%skek202412_22p5_gap_4V_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str()));
        TFile* file_blk_7v = TFile::Open(Form("%skek202412_22p5_gap_7V_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str()));
        TFile* file_blk_10v = TFile::Open(Form("%skek202412_22p5_gap_10V_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str()));

        TH1D* h_blk_4v = (TH1D*)file_blk_4v->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
        TH1D* h_blk_7v = (TH1D*)file_blk_7v->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
        TH1D* h_blk_10v = (TH1D*)file_blk_10v->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");

        plot_ExperimentData::set_1DStyle(h_blk_4v);
        plot_ExperimentData::set_1DStyle(h_blk_7v);
        plot_ExperimentData::set_1DStyle(h_blk_10v);
        h_blk_4v->Rebin(10);
        h_blk_7v->Rebin(10);
        h_blk_10v->Rebin(10);
        h_blk_4v->Scale(1/h_blk_4v->GetMaximum());
        h_blk_7v->Scale(1/h_blk_7v->GetMaximum());
        h_blk_10v->Scale(1/h_blk_10v->GetMaximum());

        h_blk_4v->SetMarkerColor(kCyan+2);
        h_blk_4v->SetLineColor(kCyan+2);
        h_blk_7v->SetMarkerColor(kBlue);
        h_blk_7v->SetLineColor(kBlue);
        h_blk_10v->SetMarkerColor(kBlue+2);
        h_blk_10v->SetLineColor(kBlue+2);

        h_blk_4v->Draw("SAMEPE");
        h_blk_7v->Draw("SAMEPE");
        h_blk_10v->Draw("SAMEPE");
        h_blk_4v->Draw("SAMEC");
        h_blk_7v->Draw("SAMEC");
        h_blk_10v->Draw("SAMEC");

        TLegend* legend_check = new TLegend(0.65, 0.50, 0.85, 0.78);
        legend_check->SetFillStyle(0);
        legend_check->SetTextSize(0.04);
        legend_check->SetBorderSize(0);

        legend_check->AddEntry(h_blk_4v, "GAP 4V", "PE");
        legend_check->AddEntry(h_blk_7v, "GAP 7V", "PE");
        legend_check->AddEntry(h_blk_10v, "GAP 10V", "PE");
        legend_check->AddEntry(h_std, "STD 10V", "PE");
        legend_check->AddEntry(h_gap, "BLK 10V", "PE");

        title.DrawLatexNDC(0.40, 0.89, "Cluster Charge");
        condition.DrawLatexNDC(0.40, 0.85, "Electron 3GeV/c @KEK");
        condition.DrawLatexNDC(0.40, 0.82, "p22p5/SeedThd60/NeighborThd60");
        condition.DrawLatexNDC(0.40, 0.79, Form("Plotted on %s", TIME_.c_str()));
        legend_check->Draw();

        canvas->SaveAs("plot/clusterCharge_blk_check.root");
    }

    bool cluster_charge_plots = true;
    if(cluster_charge_plots) {
        std::vector<std::string> chip_types = {"std", "blk", "gap"};
        std::vector<TH1D*> hists_to_delete;

        for(int i = 0; i < chip_types.size(); i++) {
            canvas->Clear();
            canvas->SetTopMargin(0.062);
            canvas->SetBottomMargin(0.14);
            canvas->SetLeftMargin(0.13);
            canvas->SetRightMargin(0.07);

            // --- 凡例の作成 ---
            // メインプロットとノイズのみ凡例に追加します。
            TLegend* legend = new TLegend(0.60, 0.65, 0.9, 0.90);
            //legend->SetHeaderText(Form("Chip: %s", CHIP_TYPE_[i].c_str()));
            legend->SetFillStyle(0); // 背景を透明に
            legend->SetBorderSize(0); // 枠線をなしに
            legend->SetTextSize(0.05);

            TLegend* legend_sub = new TLegend(0.60, 0.45, 0.9, 0.65);
            legend_sub->SetFillStyle(0);
            legend_sub->SetBorderSize(0);
            legend_sub->SetTextSize(0.05);

            //  --- ノイズのプロット ---
            TFile* noiseFile = TFile::Open(Form("%skek202412_22p5_%s_10V_SeedThd60e_NeighborThd60e_noise.root", data_dir_path.c_str(), chip_types[i].c_str()));
            if (noiseFile && !noiseFile->IsZombie()) {
                TH1D* hNoise = (TH1D*)noiseFile->Get("ClusteringAnalog/CE65_3/clusterCharge");
                if (hNoise) {
                    TH1D* hNoiseClone = (TH1D*)hNoise->Clone(Form("hNoiseClone_%s", chip_types[i].c_str()));
                    hists_to_delete.push_back(hNoiseClone); // 解放リストに追加

                    hNoiseClone->Rebin(10);
                    hNoiseClone->Scale(1.0 / hNoiseClone->GetMaximum());
                    hNoiseClone->SetStats(0);
                    hNoiseClone->SetTitle(";charge [ADC];normalized counts");
                    hNoiseClone->GetXaxis()->SetRangeUser(60, 4000);
                    hNoiseClone->GetYaxis()->SetRangeUser(0, 1.2); // 上限に少し余裕を持たせる
                    hNoiseClone->GetXaxis()->SetTitleSize(0.05); hNoiseClone->GetXaxis()->SetLabelSize(0.04); hNoiseClone->GetXaxis()->SetTitleOffset(0.8);
                    hNoiseClone->GetYaxis()->SetTitleSize(0.05); hNoiseClone->GetYaxis()->SetLabelSize(0.04); hNoiseClone->GetYaxis()->SetTitleOffset(0.8);

                    hNoiseClone->SetMarkerColor(kGray+1);
                    hNoiseClone->SetLineColor(kGray+1);
                    hNoiseClone->SetMarkerStyle(20+i);
                    hNoiseClone->SetMarkerSize(0.8);
                    //hNoiseClone->SetLineStyle(3);
                    hNoiseClone->SetLineWidth(1);
                    hNoiseClone->Draw("PE"); // 最初に描画

                    // TF1* fNoise = plot_histogram::optimise_hist_gaus(hNoiseClone, kGray);
                    // fNoise->SetLineWidth(2);
                    TF1* fNoise = new TF1(Form("fNoise_%s", chip_types[i].c_str()), "[0]*exp([1]*x)", 60, 500);
                    // fNoise->SetParameter(0, 1);
                    // fNoise->SetParameter(1, -0.01);
                    hNoiseClone->Fit(fNoise, "NLS+", "", 60, 600);
                    fNoise->SetLineColor(kGray);
                    fNoise->SetLineWidth(2);
                    fNoise->SetLineStyle(3);
                    fNoise->SetRange(hNoiseClone->GetXaxis()->GetXmin(), hNoiseClone->GetXaxis()->GetXmax());
                    fNoise->Draw("SAME");

                    legend_sub->AddEntry(hNoiseClone, Form("%s, pedestal", chip_types[i].c_str()), "pel");
                }
                //noiseFile->Close();
            }

            // --- サブチップのプロット（背景として先に描画） ---
            for (int j = 0; j < chip_types.size(); j++) {
                if (i == j) continue; // メインチップ自身はスキップ
                TFile* subFile = TFile::Open(Form("%skek202412_22p5_%s_10V_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str(), chip_types[j].c_str()));
                if (subFile && !subFile->IsZombie()) {
                    TH1D* hSub = (TH1D*)subFile->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
                    if (hSub) {
                        TH1D* hSubClone = (TH1D*)hSub->Clone(Form("hSubClone_%s_10v", chip_types[j].c_str()));
                        hists_to_delete.push_back(hSubClone); // 解放リストに追加

                        hSubClone->Rebin(10);
                        if (hSubClone->GetMaximum() > 0) hSubClone->Scale(1.0 / hSubClone->GetMaximum());

                        // サブチップは薄い色で統一
                        //int color = kSpring - 3 - j;
                        hSubClone->SetMarkerColor(kGray+2);
                        hSubClone->SetLineColor(kGray+2);
                        hSubClone->SetMarkerStyle(24+j);
                        hSubClone->SetMarkerSize(0.8);
                        hSubClone->SetLineWidth(1);

                        TF1* fSub = plot_histogram::optimise_hist_langaus(hSubClone, kGray);
                        fSub->SetLineStyle(2);
                        fSub->SetLineWidth(2);

                        hSubClone->Draw("SAME PE"); // 線なしでマーカーのみ
                        legend_sub->AddEntry(hSubClone, Form("%s, 10 V", chip_types[j].c_str()), "pel");
                    }
                        //subFile->Close();
                }
                // }

                // for (int k = 0; k < VOLTAGE_.size(); k++) {
                //     TFile* subFile = TFile::Open(Form("%skek202412_22p5_%s_%sV_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str(), chip_types[j].c_str(), VOLTAGE_[k].c_str()));
                //     if (subFile && !subFile->IsZombie()) {
                //         TH1D* hSub = (TH1D*)subFile->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
                //         if (hSub) {
                //             TH1D* hSubClone = (TH1D*)hSub->Clone(Form("hSubClone_%s_%s", chip_types[j].c_str(), VOLTAGE_[k].c_str()));
                //             hists_to_delete.push_back(hSubClone); // 解放リストに追加

                //             hSubClone->Rebin(10);
                //             if (hSubClone->GetMaximum() > 0) hSubClone->Scale(1.0 / hSubClone->GetMaximum());

                //             // サブチップは薄い色で統一
                //             int color = kSpring - 3 - k;
                //             hSubClone->SetMarkerColor(color);
                //             hSubClone->SetLineColor(color);
                //             hSubClone->SetMarkerStyle(20);
                //             hSubClone->SetMarkerSize(0.8);
                //             hSubClone->SetLineWidth(1);

                //             hSubClone->Draw("SAME PE"); // 線なしでマーカーのみ
                //         }
                //         //subFile->Close();
                //     }
                // }
            }

            // --- メインチップのプロット（一番手前に描画） ---
            for (int j = 0; j < VOLTAGE_.size(); j++) {
                TFile* mainFile = TFile::Open(Form("%skek202412_22p5_%s_%sV_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str(), chip_types[i].c_str(), VOLTAGE_[j].c_str()));
                if (mainFile && !mainFile->IsZombie()) {
                    TH1D* hMain = (TH1D*)mainFile->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");
                    if (hMain) {
                        TH1D* hMainClone = (TH1D*)hMain->Clone(Form("hMainClone_%s_%s", chip_types[i].c_str(), VOLTAGE_[j].c_str()));
                        hists_to_delete.push_back(hMainClone); // 解放リストに追加

                        hMainClone->Rebin(10);
                        if (hMainClone->GetMaximum() > 0) hMainClone->Scale(1.0 / hMainClone->GetMaximum());

                        int color = kAzure - 3 + (2 * j); // 色のバリエーション
                        hMainClone->GetXaxis()->SetRangeUser(60, 4000);
                        hMainClone->SetMarkerColor(color);
                        hMainClone->SetLineColor(color);
                        hMainClone->SetMarkerStyle(20+i); // メインとサブでマーカーを変える
                        hMainClone->SetMarkerSize(1.2);
                        if(i == 2) {
                            hMainClone->SetMarkerSize(1.4);
                        }
                        hMainClone->SetLineStyle(1);
                        hMainClone->SetLineWidth(2);
                        hMainClone->Draw("SAME PE"); // マーカー、線、エラーバーを描画

                        // TF1* fMain = new TF1(Form("fMain_%s_%s", chip_types[i].c_str(), VOLTAGE_[j].c_str()), langaufun, 60, 5000);
                        // fMain->SetLineWidth(2);
                        // fMain->SetLineColor(color);
                        // fit->SetParameters()
                        TF1* fMain = plot_histogram::optimise_hist_langaus(hMainClone, color);
                        //TF1* fMain = plot_histogram::optimise_hist_langau_expo(hMainClone, color);
                        //fMain->SetLineColorAlpha(color, 1);
                        fMain->SetLineStyle(1);
                        fMain->SetLineWidth(2);

                        //hMainClone->Draw("SAME PE"); // マーカー、線、エラーバーを描画
                        legend->AddEntry(hMainClone, Form("%s, %s V", chip_types[i].c_str(), VOLTAGE_[j].c_str()), "epl");
                    }
                    //mainFile->Close();
                }
            }

            // --- 凡例とプロットの保存 ---
            legend->Draw();
            legend_sub->Draw();
            condition.DrawLatexNDC(0.15, 0.90, "e^{-} 3GeV/c @KEK (Dec. 2024)");
            condition.DrawLatexNDC(0.15, 0.87, Form("Plotted on %s", TIME_.c_str()));
            canvas->SaveAs(Form("plot/ClusterCharge_for_%s_with_Others.pdf", chip_types[i].c_str()));
        } // end of the main loop
    } // if(cluster_charge_plots)

    // bool cluster_charge_plots = true;
    // if(cluster_charge_plots) {
    //     TFile* tmpMainROOTFile = nullptr;
    //     TFile* tmpSubROOTFile = nullptr;
    //     TFile* tmpNoiseROOTFile = nullptr;
    //     TH1D* tmpMainClusterCharge = nullptr;
    //     TH1D* tmpSubClusterCharge = nullptr;
    //     TH1D* tmpNoiseCharge = nullptr;

    //     std::vector<TH1D*> vMainClusterCharge = {};
    //     std::vector<std::vector<TH1D*>> vSubClusterCharge = {};        
    //     for(int i=0; i<CHIP_TYPE_.size(); i++) {
    //         canvas->Clear();
    //         canvas->SetTopMargin(0.062);
    //         canvas->SetBottomMargin(0.14);
    //         canvas->SetLeftMargin(0.13);
    //         canvas->SetRightMargin(0.07);
    //         TLegend* legend_main = new TLegend(0.5,0.5,0.7,0.7);

    //         tmpNoiseROOTFile = TFile::Open(Form("%skek202412_22p5_%s_10V_SeedThd60e_NeighborThd60e_noise.root", data_dir_path.c_str(), CHIP_TYPE_[i].c_str()));
    //         tmpNoiseCharge = (TH1D*)tmpNoiseROOTFile->Get("ClusteringAnalog/CE65_3/clusterCharge");
    //         //std::cout << "Mean : " << tmpNoiseCharge->GetMean() << std::endl;

    //         tmpNoiseCharge->Rebin(10);
    //         tmpNoiseCharge->Scale(1.0/tmpNoiseCharge->GetMaximum());
    //         tmpNoiseCharge->SetStats(0);
    //         tmpNoiseCharge->SetTitle(";charge [ADC];#counts");
    //         tmpNoiseCharge->GetXaxis()->SetLimits(60, 8000);
    //         tmpNoiseCharge->GetYaxis()->SetLimits(0, 1.1);
    //         tmpNoiseCharge->GetXaxis()->SetTitleSize(0.06); tmpNoiseCharge->GetXaxis()->SetLabelSize(0.04); tmpNoiseCharge->GetYaxis()->SetTitleOffset(0.7);
    //         tmpNoiseCharge->GetYaxis()->SetTitleSize(0.06); tmpNoiseCharge->GetYaxis()->SetLabelSize(0.04); tmpNoiseCharge->GetYaxis()->SetTitleOffset(0.7);

    //         tmpNoiseCharge->SetMarkerColor(kGray);
    //         tmpNoiseCharge->SetLineColor(kGray);
    //         tmpNoiseCharge->SetMarkerStyle(20);
    //         tmpNoiseCharge->SetMarkerSize(0.8);
    //         tmpNoiseCharge->SetLineStyle(1);
    //         tmpNoiseCharge->SetLineWidth(2);

    //         tmpNoiseCharge->Draw("PE");
    //         tmpNoiseCharge->Draw("C");

    //         // cluster charge of the main chip
    //         for(int j=0; j<VOLTAGE_.size(); j++) {
    //             tmpMainROOTFile = TFile::Open(Form("%skek202412_22p5_%s_%sV_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str(), CHIP_TYPE_[i].c_str(), VOLTAGE_[j].c_str()));
    //             tmpMainClusterCharge = (TH1D*)tmpMainROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");

    //             // Setting for histogram style
    //             tmpMainClusterCharge->Rebin(10);
    //             tmpMainClusterCharge->Scale(1/tmpMainClusterCharge->GetMaximum());

    //             int tmp_color = 2*j;
    //             tmpMainClusterCharge->SetMarkerColor(kAzure+tmp_color);
    //             tmpMainClusterCharge->SetMarkerColor(kAzure+tmp_color);
    //             tmpMainClusterCharge->SetMarkerStyle(20);
    //             tmpMainClusterCharge->SetMarkerSize(1.1);
    //             tmpMainClusterCharge->SetLineStyle(1);
    //             tmpMainClusterCharge->SetLineWidth(2);

    //             tmpMainClusterCharge->Draw("SAME PEC");
    //         } // VOLTAGE_

    //         // cluster charge of the sub chips
    //         for(int j=0; j<CHIP_TYPE_.size(); j++) {
    //             if(CHIP_TYPE_[i] == CHIP_TYPE_[j]) {
    //                     continue;
    //                 } else {
    //                     for(int k=0; k<VOLTAGE_.size(); k++) {
    //                         tmpSubROOTFile = TFile::Open(Form("%skek202412_22p5_%s_%sV_SeedThd60e_NeighborThd60e.root", data_dir_path.c_str(), CHIP_TYPE_[i].c_str(), VOLTAGE_[j].c_str()));
    //                         tmpSubClusterCharge = (TH1D*)tmpSubROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterCharge");

    //                         tmpSubClusterCharge->Rebin(10);
    //                         tmpSubClusterCharge->Scale(1/tmpSubClusterCharge->GetMaximum());

    //                         tmpSubClusterCharge->SetMarkerColor(kRed + j);
    //                         tmpSubClusterCharge->SetLineColor(kRed+j);
    //                         tmpSubClusterCharge->SetMarkerStyle(20);
    //                         tmpSubClusterCharge->SetMarkerSize(1.1);
    //                         tmpSubClusterCharge->SetLineStyle(1);
    //                         tmpSubClusterCharge->SetLineWidth(2);

    //                         tmpSubClusterCharge->Draw("SAME PEC");
    //                     } // VOLTAGE_
    //                 } // if chip_type
    //         } // sub CHIP_TYPE_

    //         canvas->SaveAs(Form("plot/ClusterCharge_for_%s_with_Noise_and_other.pdf",CHIP_TYPE_[i].c_str()));
    //     } // main CHIP_TYPE_
    // } // if cluster charge plots

    bool blk_reso_clsize = true;
    if(blk_reso_clsize) {
        TFile* tmpROOTFile = nullptr;
        TH1D* tmpClusterSizeTH1D = nullptr;
        TH1D* tmpResidualTH1D = nullptr;
        TF1* tmpResidualTF1 = nullptr;

        // blk
        std::vector<std::string> blk_neighbor_thresholds_1 = {"50", "60", "70", "80", "90", "100", "125", "150", "200", "250", "300"};
        std::vector<std::string> blk_neighbor_thresholds_2 = {};
        std::vector<std::string> blk_all_thresholds = blk_neighbor_thresholds_1;
        blk_all_thresholds.insert(blk_all_thresholds.end(), blk_neighbor_thresholds_2.begin(), blk_neighbor_thresholds_2.end());
        std::vector<double> blk_thresholds = {};
        std::vector<double> blk_resolution = {};
        std::vector<double> blk_clusterSize = {};
        for(int i=0; i<blk_neighbor_thresholds_1.size(); i++) {
            tmpROOTFile = TFile::Open(Form("%skek202412_22p5_blk_10V_SeedThd300e_NeighborThd%se.root", data_dir_path.c_str(), blk_neighbor_thresholds_1[i].c_str()));
            tmpClusterSizeTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterSize");
            tmpResidualTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/local_residuals/residualsX");
            
            tmpResidualTF1 = new TF1("tmpResidualTF1","gaus",-50,50);
            tmpResidualTH1D->Fit(tmpResidualTF1, "RQ");

            blk_clusterSize.push_back(tmpClusterSizeTH1D->GetMean());
            blk_resolution.push_back(tmpResidualTF1->GetParameter(2));

            blk_thresholds.push_back(std::stod(blk_neighbor_thresholds_1[i]));
        }

        // for(int i=0; i<blk_neighbor_thresholds_2.size(); i++) {
        //     tmpROOTFile = TFile::Open(Form("%skek202412_22p5_blk_10V_SeedThd%se_NeighborThd%se.root", data_dir_path.c_str(), blk_neighbor_thresholds_2[i].c_str(), blk_neighbor_thresholds_2[i].c_str()));
        //     tmpClusterSizeTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterSize");
        //     tmpResidualTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/local_residuals/residualsX");
            
        //     tmpResidualTF1 = new TF1("tmpResidualTF1","gaus",-50,50);
        //     tmpResidualTF1->Fit(tmpResidualTH1D, "RQ");

        //     blk_clusterSize.push_back(tmpClusterSizeTH1D->GetMean());
        //     blk_resolution.push_back(tmpResidualTF1->GetParameter(2));

        //     blk_thresholds.push_back(std::stod(blk_neighbor_thresholds_2[i]));
        // }


        std::vector<std::string> neighbor_thresholds = {"50", "60", "70", "80", "90", "100", "120", "150", "200", "250", "300", "350", "400", "450", "500"};
        std::vector<double> thresholds = {};
        // std
        std::vector<double> std_resolution = {};
        std::vector<double> std_clusterSize = {};
        for(int i=0; i<neighbor_thresholds.size(); i++) {
            tmpROOTFile = TFile::Open(Form("%skek202412_22p5_std_10V_SeedThd500e_NeighborThd%se.root", data_dir_path.c_str(), neighbor_thresholds[i].c_str()));
            tmpClusterSizeTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterSize");
            tmpResidualTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/local_residuals/residualsX");
            
            tmpResidualTF1 = new TF1("tmpResidualTF1","gaus",-50,50);
            tmpResidualTH1D->Fit(tmpResidualTF1, "RQ");

            std_clusterSize.push_back(tmpClusterSizeTH1D->GetMean());
            std_resolution.push_back(tmpResidualTF1->GetParameter(2));

            thresholds.push_back(std::stod(neighbor_thresholds[i]));
        }

        // gap
        std::vector<double> gap_resolution = {};
        std::vector<double> gap_clusterSize = {};
        for(int i=0; i<neighbor_thresholds.size(); i++) {
            tmpROOTFile = TFile::Open(Form("%skek202412_22p5_gap_10V_SeedThd500e_NeighborThd%se.root", data_dir_path.c_str(), neighbor_thresholds[i].c_str()));
            tmpClusterSizeTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/cluster/clusterSize");
            tmpResidualTH1D = (TH1D*)tmpROOTFile->Get("AnalysisCE65/CE65_3/local_residuals/residualsX");
            
            tmpResidualTF1 = new TF1("tmpResidualTF1","gaus",-50,50);
            tmpResidualTH1D->Fit(tmpResidualTF1, "RQ");

            gap_clusterSize.push_back(tmpClusterSizeTH1D->GetMean());
            gap_resolution.push_back(tmpResidualTF1->GetParameter(2));
        }

        TGraph* gr_res_blk = new TGraph(blk_thresholds.size(), &blk_thresholds[0], &blk_resolution[0]);
        TGraph* gr_cs_blk  = new TGraph(blk_thresholds.size(), &blk_thresholds[0], &blk_clusterSize[0]);
        TGraph* gr_res_std = new TGraph(thresholds.size(), &thresholds[0], &std_resolution[0]);
        TGraph* gr_cs_std  = new TGraph(thresholds.size(), &thresholds[0], &std_clusterSize[0]);
        TGraph* gr_res_gap = new TGraph(thresholds.size(), &thresholds[0], &gap_resolution[0]);
        TGraph* gr_cs_gap  = new TGraph(thresholds.size(), &thresholds[0], &gap_clusterSize[0]);

        // Style Resolution graphs (blue tones, different markers)
        gr_res_blk->SetMarkerStyle(20); gr_res_blk->SetMarkerColorAlpha(kAzure+2,   1); gr_res_blk->SetLineColorAlpha(kAzure+2,   1); gr_res_blk->SetLineStyle(1); gr_res_blk->SetLineWidth(2); gr_res_blk->SetMarkerSize(1.1);
        gr_res_std->SetMarkerStyle(21); gr_res_std->SetMarkerColorAlpha(kAzure-4,   1); gr_res_std->SetLineColorAlpha(kAzure-4, 0.6); gr_res_std->SetLineStyle(2); gr_res_std->SetLineWidth(2); gr_res_std->SetMarkerSize(1.1);
        gr_res_gap->SetMarkerStyle(22); gr_res_gap->SetMarkerColorAlpha(kAzure-7, 0.6); gr_res_gap->SetLineColorAlpha(kAzure-7, 0.6); gr_res_gap->SetLineStyle(3); gr_res_gap->SetLineWidth(2); gr_res_gap->SetMarkerSize(1.3);
        
        // Style Cluster Size graphs (red tones, different markers)
        gr_cs_blk->SetMarkerStyle(20); gr_cs_blk->SetMarkerColorAlpha(kRed+1,   1); gr_cs_blk->SetLineColorAlpha(kRed+1,   1); gr_cs_blk->SetLineStyle(1); gr_cs_blk->SetLineWidth(2); gr_cs_blk->SetMarkerSize(1.1);
        gr_cs_std->SetMarkerStyle(21); gr_cs_std->SetMarkerColorAlpha(kRed-4, 0.7); gr_cs_std->SetLineColorAlpha(kRed-4, 0.6); gr_cs_std->SetLineStyle(2); gr_cs_std->SetLineWidth(2); gr_cs_std->SetMarkerSize(1.1);
        gr_cs_gap->SetMarkerStyle(22); gr_cs_gap->SetMarkerColorAlpha(kRed-7,   1); gr_cs_gap->SetLineColorAlpha(kRed-7, 0.6); gr_cs_gap->SetLineStyle(3); gr_cs_gap->SetLineWidth(2); gr_cs_gap->SetMarkerSize(1.3);

        // --- 3. Plotting ---
        canvas->Clear();
        canvas->SetTopMargin(0.052);
        canvas->SetBottomMargin(0.145);
        canvas->SetLeftMargin(0.13);
        canvas->SetRightMargin(0.13);

        // Pad for left axis (Resolution)
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
        pad1->SetGrid();
        pad1->Draw();
        pad1->cd();

        TMultiGraph* mg_resolution = new TMultiGraph();
        //mg_resolution->Add(gr_res_blk, "PL");
        mg_resolution->Add(gr_res_std, "PL");
        mg_resolution->Add(gr_res_gap, "PL");
        mg_resolution->Add(gr_res_blk, "PL");
        
        mg_resolution->Draw("A");
        mg_resolution->SetTitle(";threshold [ADC];resolution in x [um]");
        mg_resolution->GetXaxis()->SetLimits(40, 510); // Set X-axis range
        mg_resolution->GetYaxis()->SetRangeUser(3.5, 9.5); // Set Y-axis range for resolution

        // Axis styling
        mg_resolution->GetXaxis()->SetTitleSize(0.06);
        mg_resolution->GetXaxis()->SetLabelSize(0.04);
        mg_resolution->GetXaxis()->SetTitleOffset(0.6);
        mg_resolution->GetYaxis()->SetTitleSize(0.06);
        mg_resolution->GetYaxis()->SetLabelSize(0.04);
        mg_resolution->GetYaxis()->SetTitleOffset(0.7);
        //mg_resolution->GetYaxis()->SetTitleColor(kAzure-4);
        //mg_resolution->GetYaxis()->SetLabelColor(kAzure-4);

        // Transparent pad for right axis (Cluster Size)
        canvas->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
        pad2->SetFillStyle(0); // Transparent
        pad2->SetFrameFillStyle(0);
        pad2->Draw();
        pad2->cd();

        TMultiGraph* mg_clustersize = new TMultiGraph();
        //mg_clustersize->Add(gr_cs_blk, "PL");
        mg_clustersize->Add(gr_cs_std, "PL");
        mg_clustersize->Add(gr_cs_gap, "PL");
        mg_clustersize->Add(gr_cs_blk, "PL");
        
        // Draw with invisible axes to align with pad1
        mg_clustersize->Draw("A"); 
        mg_clustersize->GetXaxis()->SetLimits(40, 510);
        mg_clustersize->GetYaxis()->SetRangeUser(0.5, 6.5);
        mg_clustersize->GetXaxis()->SetLabelSize(0);
        mg_clustersize->GetXaxis()->SetTickLength(0);
        mg_clustersize->GetYaxis()->SetLabelSize(0);
        mg_clustersize->GetYaxis()->SetTickLength(0);


        // Right Y-axis (TGaxis)
        pad2->Update(); // Important: update pad before getting coordinates
        double y_min_cs = 0.5;
        double y_max_cs = 6.5;
        TGaxis *axis_cs = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), y_min_cs, y_max_cs, 510, "+L");
        axis_cs->SetTitle("mean cluster size");
        //axis_cs->SetTitleColor(kRed+1);
        //axis_cs->SetLineColor(kRed+1);
        //axis_cs->SetLabelColor(kRed+1);
        axis_cs->SetTitleSize(0.06);
        axis_cs->SetLabelSize(0.04);
        axis_cs->SetTitleOffset(0.7);
        axis_cs->SetLabelFont(42);
        axis_cs->SetTitleFont(42);
        axis_cs->Draw();

        // Legend
        canvas->cd();
        TLegend* legend_comp_residual = new TLegend(0.70, 0.4, 0.93, 0.63);
        legend_comp_residual->SetFillStyle(0);
        legend_comp_residual->SetBorderSize(0);
        legend_comp_residual->SetTextSize(0.04);
        //legend_comp_residual->SetHeader("resolution", "C");
        legend_comp_residual->AddEntry(gr_res_blk, "blk, 10V", "pl");
        legend_comp_residual->AddEntry(gr_res_std, "std, 10V", "pl");
        legend_comp_residual->AddEntry(gr_res_gap, "gap, 10V", "pl");
        legend_comp_residual->Draw();

        TLegend* legend_comp_clsize = new TLegend(0.65, 0.4, 0.88, 0.63);
        legend_comp_clsize->SetFillStyle(0);
        legend_comp_clsize->SetBorderSize(0);
        legend_comp_clsize->AddEntry(gr_cs_blk, " ", "pl");
        legend_comp_clsize->AddEntry(gr_cs_std, " ", "pl");
        legend_comp_clsize->AddEntry(gr_cs_gap, " ", "pl");
        legend_comp_clsize->Draw();
        

        condition.DrawLatexNDC(0.12, 0.87, "e^{-} 3GeV/c @KEK (Dec. 2024)");
        condition.DrawLatexNDC(0.65, 0.87, Form("Plotted on %s", TIME_.c_str()));

        // Save the plot
        canvas->SaveAs("plot/Resolution_vs_ClusterSize_Comparison_forBLK.pdf");
    }

    bool seed_charge = true;
    if(seed_charge) {
        canvas->Clear();
    }


    output->Close();
}

void plot_ExperimentData::run_inPixel() {
    LOG_STATUS.source("plot_ExperimentData::run_inPixel") << "Start run for in-pixel analysis.";
    
    std::string output_file_name = "/home/towa/alice3/plotter/plot/experimentalData_inPixel.root";
    std::string data_dir_path = "/home/towa/alice3/hist/kek202412/";
    TFile* output = TFile::Open(output_file_name.c_str(), "RECREATE");

    TDirectory* clsize = output->mkdir("cluster_size");
    TDirectory* clusterCharge = output->mkdir("cluster_charge");
    TDirectory* seedCharge = output->mkdir("seed_charge");
    TDirectory* residual = output->mkdir("residual");

    std::vector<TDirectory*> outputDir = {clsize, clusterCharge, seedCharge, residual};

    TFile* inputROOTFile;
    TProfile2D* p_clSize = nullptr;
    TProfile2D* p_clusterCharge = nullptr;
    TProfile2D* p_seedCharge = nullptr;
    TProfile2D* p_residual = nullptr;
    TH2D* hClSize = nullptr;
    TH2D* hClusterCharge = nullptr;
    TH2D* hSeedCharge = nullptr;
    TH2D* hResidual = nullptr;

    TCanvas* canvas = new TCanvas("canvas","canvas",800,600);
    TLatex title;
    TLatex condition;
    title.SetTextSize(0.04);
    title.SetTextFont(62);
    condition.SetTextSize(0.03);
    condition.SetTextFont(62);

    gStyle->SetPalette(55);

    std::string chip_variation;
    std::string chip_variation_text;

    for(int i=0; i<PIXEL_PITCH_.size(); i++) {
        for(int j=0; j<CHIP_TYPE_.size(); j++) {
            for(int k=0; k<VOLTAGE_.size(); k++) {
                for(int l=0; l<SEED_THRESHOLD_.size(); l++) {
                    for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
                        inputROOTFile = TFile::Open(Form("%skek202412_%s_%s_%sV_SeedThd%se_NeighborThd%se.root", data_dir_path.c_str(), PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str()));
                        chip_variation = Form("%s_%s_%sV_SeedThd%se_NeighborThd%se", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str());
                        chip_variation_text = Form("p%s/%s/%sV/SeedThd%s/NeighborThd%s", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str(), NEIGHBOR_THRESHOLD_[n].c_str());
                        gStyle->SetPalette(kViridis);

                        p_clSize = (TProfile2D*)inputROOTFile->Get("AnalysisCE65/CE65_3/npxvsxmym");
                        p_clusterCharge = (TProfile2D*)inputROOTFile->Get("AnalysisCE65/CE65_3/qvsxmym");
                        p_seedCharge = (TProfile2D*)inputROOTFile->Get("AnalysisCE65/CE65_3/pxqvsxmym");
                        p_residual = (TProfile2D*)inputROOTFile->Get("AnalysisCE65/CE65_3/rmsxyvsxmym");

                        hClSize = convert_toTH2D(p_clSize);
                        hClusterCharge = convert_toTH2D(p_clusterCharge);
                        hSeedCharge = convert_toTH2D(p_seedCharge);
                        hResidual = convert_toTH2D(p_residual);
                    
                        hClSize->SetMinimum(0);
                        hClSize->SetMaximum(8);
                        hClSize->SetTitle(";x w/in pixel [um];y w/in pixel [um];cluster size");
                        plot_ExperimentData::set_2DSURFStyle(canvas, hClSize);
                        title.DrawLatexNDC(0.20, 0.89, "Mean of the Cluster Size");
                        condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                        condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                        condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                        clsize->cd();
                        canvas->Write(Form("clsize_ce65_%s", chip_variation.c_str()));
                        //canvas->SaveAs(Form("inPixel/experimentalData_clsize_ce65_%s.pdf", chip_variation.c_str()));
                    
                        hResidual->Scale(1000.0);
                        hResidual->SetMinimum(0);
                        hResidual->SetMaximum(15);
                        hResidual->SetTitle(";x w/in pixel [um];y w/in pixel [um];r_{exp}-r_{hit} [um]");
                        plot_ExperimentData::set_2DSURFStyle(canvas, hResidual);
                        title.DrawLatexNDC(0.20, 0.89, "Mean of the Residual (r_{exp} - r_{hit})");
                        condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                        condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                        condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                        residual->cd();
                        canvas->Write(Form("residual_ce65_%s", chip_variation.c_str()));

                        hClusterCharge->SetMinimum(0);
                        hClusterCharge->SetMaximum(3000);
                        hClusterCharge->SetTitle(";x w/in pixel [um];y w/in pixel [um];cluster charge [adu]");
                        plot_ExperimentData::set_2DSURFStyle(canvas, hClusterCharge);
                        title.DrawLatexNDC(0.20, 0.89, "Mean of the Cluster Charge");
                        condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                        condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                        condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                        clusterCharge->cd();
                        canvas->Write(Form("clusterCharge_ce65_%s", chip_variation.c_str()));

                        hSeedCharge->SetMinimum(0);
                        hSeedCharge->SetMaximum(1600);
                        hSeedCharge->SetTitle(";x w/in pixel [um];y w/in pixel [um];seed charge [adu]");
                        plot_ExperimentData::set_2DSURFStyle(canvas, hSeedCharge);
                        title.DrawLatexNDC(0.20, 0.89, "Mean of the Seed Charge");
                        condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                        condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                        condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                        seedCharge->cd();
                        canvas->Write(Form("seedCharge_ce65_%s", chip_variation.c_str()));

                        inputROOTFile->Close();
                    } // NEIGHBOR_THRESHOLD_
                } // SEED_THRESHOLD_
            } // VOLTAE_
        } // CHIP_TYPE_
    } // PIXEL_PITCH_

    output->Close();

    plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), "cluster_size", Form("plot/clusterSize_ce65_inPixel.pdf"));
    plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), "residual", Form("plot/residual_ce65_inPixel.pdf"));
    plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), "cluster_charge", Form("plot/clusterCharge_ce65_inPixel.pdf"));
    plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), "seed_charge", Form("plot/seedCharge_ce65_inPixel.pdf"));

}