#include "plot_ExperimentData.h"

plot_ExperimentData::plot_ExperimentData() {
    LOG_STATUS.source("plot_ExperimentData::plot_ExperimentData") << "Plot_ExperimentData object is created.";
    PIXEL_PITCH_ = {"22p5"};
    CHIP_TYPE_ = {"std", "gap"};
    //CHIP_TYPE_ = {"blk"};
    VOLTAGE_ = {"10", "7", "4"};
    //VOLTAGE_ = {"10"};
    //SEED_THRESHOLD_ = {"400"};
    //SEED_THRESHOLD_ = {"100"};
    SEED_THRESHOLD_ = {"500"};
    NEIGHBOR_THRESHOLD_ = {"50", "60" ,"70", "80", "90", "100", "150", "200", "250", "300", "350", "400", "450", "500"}; 
    //NEIGHBOR_THRESHOLD_ = {"50","60" ,"70", "80", "90", "100"}; 
    //NEIGHBOR_THRESHOLD_ = {"60", "100", "150", "200", "300", "400"};
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
            canvas->SaveAs(Form("perPixelCharge_ce65_%s.pdf", chip_variation.c_str()));
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
                    canvas->SaveAs(Form("clusterCharge_ce65_%s.pdf", chip_variation.c_str()));

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
                    canvas->SaveAs(Form("residualX_ce65_%s.pdf", chip_variation.c_str()));
                } // SEED_THRESHOLD_
            } // VOLTAGE_
        } // CHIP_TYPE_
    } // PIXEL_PTICH_

    // for(int i=0; i<PIXEL_PITCH_.size(); i++) {
    //     for(int j=0; j<CHIP_TYPE_.size(); j++) {
    //         for(int k=0; k<VOLTAGE_.size(); k++) {
    //             for(int l=0; l<SEED_THRESHOLD_.size(); l++) {
    //                 canvas->Clear();
    //                 canvas->SetTopMargin(0.062);
    //                 canvas->SetBottomMargin(0.14);
    //                 canvas->SetLeftMargin(0.2);
    //                 canvas->SetRightMargin(0.2);

    //                 TGraph* gr_resolution = new TGraph();
    //                 TGraph* gr_clustersize = new TGraph();

    //                 plot_ExperimentData::set_GraphStyle(gr_resolution);
    //                 plot_ExperimentData::set_GraphStyle(gr_clustersize);
    //                 gr_resolution->SetMarkerStyle(20);
    //                 gr_clustersize->SetMarkerStyle(21);

    //                 for(int n=0; n<NEIGHBOR_THRESHOLD_.size(); n++) {
    //                     double x_val = std::stod(NEIGHBOR_THRESHOLD_[n]);
    //                     double y_resolution = vResolutionX[i][j][k][l][n];
    //                     double y_clustersize = vMeanClusterSize[i][j][k][l][n];

    //                     gr_resolution->SetPoint(n, x_val, y_resolution);
    //                     gr_clustersize->SetPoint(n, x_val, y_clustersize);
    //                 } // NEIGHBOR_THRESHOLD_

    //                 // Drawing
    //                 canvas->cd();
    //                 TPad* pad1 = new TPad("pad1","pad1",0,0,1,1);
    //                 pad1->Draw();
    //                 pad1->cd();
    //                 gr_resolution->Draw("APL");
    //                 gr_resolution->SetTitle(";threshold [ADC];resolution x[um]");
    //                 gr_resolution->GetXaxis()->SetTitleSize(0.06);
    //                 gr_resolution->GetXaxis()->SetLabelSize(0.04);
    //                 gr_resolution->GetXaxis()->SetTitleOffset(0.8);
    //                 gr_resolution->GetYaxis()->SetRangeUser(0, 10);

    //                 gr_resolution->GetYaxis()->SetTitleFont(42);
    //                 gr_resolution->GetYaxis()->SetLabelFont(42);
    //                 gr_resolution->GetYaxis()->SetTitleSize(0.06);
    //                 gr_resolution->GetYaxis()->SetLabelSize(0.04);
    //                 gr_resolution->GetYaxis()->SetTitleOffset(0.8);
                    
    //                 canvas->cd();
    //                 TPad* pad2 = new TPad("pad2","pad2",0,0,1,1);
    //                 pad2->SetFillStyle(0);
    //                 pad2->SetFrameFillStyle(0);
    //                 pad2->Draw();
    //                 pad2->cd();
    //                 gr_clustersize->Draw("APL");
    //                 gr_clustersize->GetYaxis()->SetRangeUser(0,6);
    //                 gr_clustersize->GetXaxis()->SetLabelSize(0);
    //                 gr_clustersize->GetXaxis()->SetTitleSize(0);
    //                 gr_clustersize->GetYaxis()->SetLabelSize(0);
    //                 gr_clustersize->GetYaxis()->SetTitleSize(0);
    //                 gr_clustersize->GetYaxis()->SetTickSize(0);                    

    //                 canvas->Update();
    //                 TGaxis *axis_clustersize = new TGaxis(
    //                 gPad->GetUxmax(), gPad->GetUymin(), 
    //                 gPad->GetUxmax(), gPad->GetUymax(),
    //                 gr_clustersize->GetYaxis()->GetXmin(), gr_clustersize->GetYaxis()->GetXmax(), 510, "+L");
    //                 axis_clustersize->SetTitle("Mean Cluster Size");
    //                 axis_clustersize->SetTitleOffset(0.8);
    //                 axis_clustersize->SetTitleSize(0.06);
    //                 axis_clustersize->SetLabelSize(0.04);
    //                 axis_clustersize->SetTitleFont(42);
    //                 axis_clustersize->SetLabelFont(42);
    //                 // axis_clustersize->SetTitleColor(kRed);
    //                 // axis_clustersize->SetLineColor(kRed);
    //                 // axis_clustersize->SetLabelColor(kRed);
    //                 //axis_clustersize->SetTitleOffset(1.2);
    //                 axis_clustersize->Draw();

    //                 canvas->cd(); 
    //                 TLegend* legend_graph = new TLegend(0.15, 0.75, 0.45, 0.88);
    //                 legend_graph->SetFillStyle(0);
    //                 legend_graph->SetBorderSize(0);
    //                 legend_graph->AddEntry(gr_resolution, "Position Resolution", "pl");
    //                 legend_graph->AddEntry(gr_clustersize, "Mean Cluster Size", "pl");
    //                 legend_graph->Draw();
                    
                    
    //                 chip_variation = Form("%s_%s_%sV_SeedThd%se", PIXEL_PITCH_[i].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[k].c_str(), SEED_THRESHOLD_[l].c_str());
    //                 canvas->Write(Form("neighborThd_Scan_%s.pdf", chip_variation.c_str()));
    //                 canvas->SaveAs(Form("NeighborThd_Scan_%s.pdf", chip_variation.c_str()));
    //             } // SEED_THRESHOLD_
    //         } // VOLTAGE_
    //     } // CHIP_TYPE_
    // } // PIXEL_PITCH_

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

        canvas->SaveAs("clusterCharge_blk_check.root");
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