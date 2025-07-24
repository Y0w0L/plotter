#include "plot_ExperimentData.h"

plot_ExperimentData::plot_ExperimentData() {
    LOG_STATUS.source("plot_ExperimentData::plot_ExperimentData") << "Plot_ExperimentData object is created.";
    PIXEL_PITCH_ = {"22p5"};
    CHIP_TYPE_ = {"std", "gap"};
    //CHIP_TYPE_ = {"std"};
    VOLTAGE_ = {"10", "7", "4"};
    //VOLTAGE_ = {"10"};
    //SEED_THRESHOLD_ = {"400"};
    //SEED_THRESHOLD_ = {"100"};
    SEED_THRESHOLD_ = {"500"};
    //NEIGHBOR_THRESHOLD_ = {"0", "10", "20", "30", "40", "50", "60" ,"70", "80", "90", "100", "125", "150", "200", "250", "300", "350", "400"}; 
    //NEIGHBOR_THRESHOLD_ = {"0", "10", "20", "30", "40", "50", "60" ,"70", "80", "90", "100"}; 
    NEIGHBOR_THRESHOLD_ = {"60", "100", "150", "200", "300", "400"};
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

void plot_ExperimentData::set_1DStyle(TH1D* hist) {
    hist->SetStats(0);
    hist->GetXaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->SetTitleFont(42, "XYZ");

    hist->SetMarkerSize(0.6);
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
                            hClusterCharge->GetXaxis()->SetRangeUser(500, 5000);
                            hClusterCharge->GetYaxis()->SetRangeUser(0, 1400);
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
                        hResidualX->SetTitle(";x_{track} - x_{hit} [um];counts");
                        hResidualX->Rebin(4);
                        plot_ExperimentData::set_1DStyle(hResidualX);
                        hResidualX->SetMarkerColorAlpha(my_colors[n], 1);
                        hResidualX->SetLineColorAlpha(my_colors[n], 1);
                        legend->AddEntry(hResidualX, Form("Thd%s", NEIGHBOR_THRESHOLD_[n].c_str()), "PE");

                        if(n==0) {
                            hResidualX->GetXaxis()->SetRangeUser(-50, 50);
                            hResidualX->GetYaxis()->SetRangeUser(0, 1000);
                            hResidualX->Draw("PE");
                        } else {
                            hResidualX->Draw("same PE");
                        }

                        fResidualX = new TF1("fResidualX", "gaus", -50, 50);
                        fResidualX->SetLineColor(my_colors[n]);
                        hResidualX->Fit(fResidualX, "RQ");
                    }// NEIGHBOR_THRESHOLD_

                    title.DrawLatexNDC(0.60, 0.89, "Cluster Charge");
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