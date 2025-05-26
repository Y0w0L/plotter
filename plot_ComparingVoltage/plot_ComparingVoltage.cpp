#include "plot_ComparingVoltage.h"

const std::vector<std::string> chip_type_ = {"std", "gap"};
const std::vector<std::string> pixel_pitch_ = {"15", "22p5"};
const std::vector<std::string> voltage_ = {"4", "7", "10"};
const std::vector<std::string> threshold_ = {"0", "10", "25", "50", "100", "200", "300", "400"};
const std::vector<std::string> model_ = {"masetti", "canali", "jacoboni"};

std::string chip_type, pixel_pitch, voltage, threshold, model;

// histogram x
const int x_Clsize_min = 1;
const int x_Clsize_max = 10;
const int y_Clsize_min = 0;
const int y_Clsize_max = 1;
const int x_Charge_min = 0;
const int x_Charge_max = 3000;
const int y_Charge_min = 0;
const int y_Charge_max = 700;
const int x_Residual_min = -15;
const int x_Residual_max = 15;
const int y_Residual_min = 0;
const int y_Residual_max = 800;

plot_ComparingVoltage::plot_ComparingVoltage() {
    LOG_STATUS.source("plot_ComparingVoltage::plot_ComparingVoltage") << "plot_ComparingVoltage object is created.";
}

void plot_ComparingVoltage::get_voltageClsize(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage) {
    const std::string hClsize_path = "DetectorHistogrammer/CE65/cluster_size/cluster_size";
    std::vector<TH1D*> vClsize = {};
    
    // Get cluster size histrogram in input_ROOTFile(vector), then append histogram to vector
    for (int i=0; i<input_ROOTFile.size(); ++i) {
        TH1D* hClsize = (TH1D*)input_ROOTFile[i]->Get(hClsize_path.c_str());
        if(hClsize == nullptr) {
            LOG_WARNING.source("plot_ComparingVoltage::get_voltageClsize") << "Histogram is nullptr";
            break;
        }
        hClsize->Scale(1/hClsize->Integral());
        vClsize.push_back(hClsize);
    }

    // Writing histogram to TCanvas
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    c->SetTopMargin(0.062);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.03);

    for(int i=0; i<vClsize.size(); ++i) {
        vClsize[i]->SetMarkerSize(1.5);
        vClsize[i]->SetMarkerStyle(20);
        if(i == 0) {
            vClsize[i]->SetTitle(";cluster size;clusters (normalized)");
            vClsize[i]->GetXaxis()->SetRangeUser(x_Clsize_min, x_Clsize_max);
            vClsize[i]->GetYaxis()->SetRangeUser(y_Clsize_min, y_Clsize_max);
            vClsize[i]->GetXaxis()->SetTitleOffset(0.8);
            vClsize[i]->GetYaxis()->SetTitleOffset(0.8);
            vClsize[i]->GetXaxis()->SetTitleSize(0.06);
            vClsize[i]->GetYaxis()->SetTitleSize(0.06);
            vClsize[i]->GetXaxis()->SetLabelSize(0.04);
            vClsize[i]->GetYaxis()->SetLabelSize(0.04);
            vClsize[i]->SetMarkerColor(1);
            vClsize[i]->Draw("histE");
        } else {
            vClsize[i]->SetMarkerColor(i+1);
            vClsize[i]->SetLineColor(i+1);
            if(i == 2) {
                vClsize[i]->SetLineColor(kGreen+2);
                vClsize[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                vClsize[i]->SetLineColor(kCyan);
                vClsize[i]->SetMarkerColor(kCyan);
            }
            vClsize[i]->Draw("samehistE");
        }
    }

    // Write TLegend on TCanvas
    const double y_legend_max = 0.72 - (0.06*vClsize.size());
    TLegend* legend = new TLegend(0.7, 0.72, 0.88, y_legend_max);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    for(int i=0; i<vClsize.size(); ++i) {
        legend->AddEntry(vClsize[i], Form("%sV", voltage[i].c_str()), "p");
    }
    legend->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatexNDC(0.63, 0.88, "Cluster size");
    latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("p%s/%s/Thd%se", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));
    if(model == "masetti") {
        latex.DrawLatexNDC(0.63, 0.73, "Masetti-Canali");
    }
    if(model == "jacoboni") {
        latex.DrawLatexNDC(0.63, 0.73, "Jacoboni-Canali");
    }
    if(model == "canali") {
        latex.DrawLatexNDC(0.63, 0.73, "Canali");
    }

    output_ROOTFile->cd();
    output_ROOTFile->cd(model.c_str());
    c->Write(Form("hClsize_%s_%s_Thd%s", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));

    // delete c;
    // delete legend;
    // delete vClsize;
    // delete y_legend_max;
    // delete hClsize_path;
}

void plot_ComparingVoltage::get_voltageCharge(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage) {
    const std::string hClusterCharge_path = "DetectorHistogrammer/CE65/charge/cluster_charge";
    const std::string hSeedCharge_path = "DetectorHistogrammer/CE65/charge/seed_charge";

    std::vector<TH1D*> vClusterCharge = {};
    std::vector<TH1D*> vSeedCharge = {};
    std::vector<TF1*> vFitClusterCharge = {};
    std::vector<TF1*> vFitSeedCharge = {};

    int nBins;
    double xMax;

    // Get ClusterCharge distribution in input_ROOTFile
    for(int i=0; i<voltage.size(); ++i) {
        TH1D* hClusterCharge = (TH1D*)input_ROOTFile[i]->Get(hClusterCharge_path.c_str());

        if(hClusterCharge == nullptr) {
            LOG_WARNING.source("plot_ComparingVoltage::get_voltageCharge") << "Histogram is nullptr.";
            break;
        }

        nBins = hClusterCharge->GetNbinsX();
        xMax = hClusterCharge->GetXaxis()->GetBinCenter(nBins)*1000;
        TH1D* hClone = (TH1D*)hClusterCharge->Clone();
        hClusterCharge = new TH1D(hClusterCharge->GetName(), ";charge [e];events", nBins, 0, xMax);
        for(int j=0; j<nBins; j++) {
            hClusterCharge->SetBinContent(j, hClone->GetBinContent(j));
        }
        vClusterCharge.push_back(hClusterCharge);
        delete hClone;
    }

    // Get SeedCharge distribution in input_ROOTFile
    for(int i=0; i<voltage.size(); ++i) {
        TH1D* hSeedCharge = (TH1D*)input_ROOTFile[i]->Get(hSeedCharge_path.c_str());
        if(hSeedCharge == nullptr) {
            LOG_WARNING.source("plot_ComparingVoltage::get_voltageCharge") << "Histogram is nullptr.";
            break;
        }

        nBins = hSeedCharge->GetNbinsX();
        xMax = hSeedCharge->GetXaxis()->GetBinCenter(nBins)*1000;
        TH1D* hClone = (TH1D*)hSeedCharge->Clone();
        hSeedCharge = new TH1D(hSeedCharge->GetName(), ";charge [e];events", nBins, 0, xMax);
        for(int j=0; j<nBins; ++j) {
            hSeedCharge->SetBinContent(j, hClone->GetBinContent(j));
        }
        vSeedCharge.push_back(hSeedCharge);
        delete hClone;
    }

    // Writing ClusterCharge plot
    // Writing histogram to TCanvas
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    c->SetTopMargin(0.062);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.03);

    for(int i=0; i<vClusterCharge.size(); ++i) {
        vClusterCharge[i]->SetMarkerSize(1.5);
        vClusterCharge[i]->SetMarkerStyle(20);
        if(i == 0) {
            vClusterCharge[i]->SetTitle(";charge [e]; events");
            vClusterCharge[i]->GetXaxis()->SetRangeUser(x_Charge_min, x_Charge_max);
            vClusterCharge[i]->GetYaxis()->SetRangeUser(y_Charge_min, y_Charge_max);
            vClusterCharge[i]->GetXaxis()->SetTitleOffset(0.8);
            vClusterCharge[i]->GetYaxis()->SetTitleOffset(0.8);
            vClusterCharge[i]->GetXaxis()->SetTitleSize(0.06);
            vClusterCharge[i]->GetYaxis()->SetTitleSize(0.06);
            vClusterCharge[i]->GetXaxis()->SetLabelSize(0.04);
            vClusterCharge[i]->GetYaxis()->SetLabelSize(0.04);
            vClusterCharge[i]->SetMarkerColor(1);
            vClusterCharge[i]->SetLineColor(1);
            vClusterCharge[i]->Draw("PE");
        } else {
            vClusterCharge[i]->SetMarkerColor(i+1);
            vClusterCharge[i]->SetLineColor(i+1);
            if(i == 2) {
                vClusterCharge[i]->SetLineColor(kGreen+2);
                vClusterCharge[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                vClusterCharge[i]->SetLineColor(kCyan);
                vClusterCharge[i]->SetMarkerColor(kCyan);
            }
            vClusterCharge[i]->Draw("samePE");
        }
    }
    vFitClusterCharge = plot_MobilityModel::get_landauFits(vClusterCharge);

    // Write TLegend on TCanvas
    const double y_legend_max = 0.72 - (0.06*vClusterCharge.size());
    const double y_fit_legend_max = 0.88 - (0.071*vClusterCharge.size());
    TLegend* legend = new TLegend(0.7, 0.72, 0.88, y_legend_max);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    for(int i=0; i<vClusterCharge.size(); ++i) {
        legend->AddEntry(vClusterCharge[i], Form("%sV", voltage[i].c_str()), "p");
    }
    legend->Draw();
    TLegend* fit_legend = new TLegend(0.12, 0.88, 0.4, y_fit_legend_max);
    fit_legend->SetFillStyle(0);
    fit_legend->SetTextSize(0.03);
    fit_legend->SetBorderSize(0);
    for(int i=0; i<vClusterCharge.size(); ++i) {
        //fit_legend->AddEntry(fits[i], names[i].c_str(), "l");
        //fit_legend->AddEntry(fits[i], Form("mean: %.2f/sigma: %.2f", fits[i]->GetParameter(1), fits[i]->GetParameter(2)), "");
    }
    fit_legend->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatexNDC(0.63, 0.88, "Cluster charge");
    latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("p%s/%s/Thd%se", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));
    if(model == "masetti") {
        latex.DrawLatexNDC(0.63, 0.73, "Masetti-Canali");
    }
    if(model == "jacoboni") {
        latex.DrawLatexNDC(0.63, 0.73, "Jacoboni-Canali");
    }
    if(model == "canali") {
        latex.DrawLatexNDC(0.63, 0.73, "Canali");
    }

    output_ROOTFile->cd();
    output_ROOTFile->cd(model.c_str());
    c->Write(Form("hClusterCharge_%s_%s_Thd%s", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));

    // Writing SeedCluster distribution in inputROOTFile to output file
    TCanvas* c_seed = new TCanvas("c", "c", 800, 600);
    c_seed->SetTopMargin(0.062);
    c_seed->SetBottomMargin(0.14);
    c_seed->SetLeftMargin(0.10);
    c_seed->SetRightMargin(0.03);

    for(int i=0; i<vClusterCharge.size(); ++i) {
        vSeedCharge[i]->SetMarkerSize(1.5);
        vSeedCharge[i]->SetMarkerStyle(20);
        if(i == 0) {
            vSeedCharge[i]->SetTitle(";charge [e]; events");
            vSeedCharge[i]->GetXaxis()->SetRangeUser(x_Charge_min, x_Charge_max);
            vSeedCharge[i]->GetYaxis()->SetRangeUser(y_Charge_min, y_Charge_max);
            vSeedCharge[i]->GetXaxis()->SetTitleOffset(0.8);
            vSeedCharge[i]->GetYaxis()->SetTitleOffset(0.8);
            vSeedCharge[i]->GetXaxis()->SetTitleSize(0.06);
            vSeedCharge[i]->GetYaxis()->SetTitleSize(0.06);
            vSeedCharge[i]->GetXaxis()->SetLabelSize(0.04);
            vSeedCharge[i]->GetYaxis()->SetLabelSize(0.04);
            vSeedCharge[i]->SetMarkerColor(1);
            vSeedCharge[i]->SetLineColor(1);
            vSeedCharge[i]->Draw("PE");
        } else {
            vSeedCharge[i]->SetMarkerColor(i+1);
            vSeedCharge[i]->SetLineColor(i+1);
            if(i == 2) {
                vSeedCharge[i]->SetLineColor(kGreen+2);
                vSeedCharge[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                vSeedCharge[i]->SetLineColor(kCyan);
                vSeedCharge[i]->SetMarkerColor(kCyan);
            }
            vSeedCharge[i]->Draw("samePE");
        }
    }
    vFitSeedCharge = plot_MobilityModel::get_landauFits(vSeedCharge);

    // Write TLegend on TCanvas
    TLegend* legend_seed = new TLegend(0.7, 0.72, 0.88, y_legend_max);
    legend_seed->SetFillStyle(0);
    legend_seed->SetTextSize(0.04);
    legend_seed->SetBorderSize(0);
    for(int i=0; i<vSeedCharge.size(); ++i) {
        legend_seed->AddEntry(vSeedCharge[i], Form("%sV", voltage[i].c_str()), "p");
    }
    legend_seed->Draw();
    TLegend* fit_legend_seed = new TLegend(0.12, 0.88, 0.4, y_fit_legend_max);
    fit_legend_seed->SetFillStyle(0);
    fit_legend_seed->SetTextSize(0.03);
    fit_legend_seed->SetBorderSize(0);
    for(int i=0; i<vSeedCharge.size(); ++i) {
        //fit_legend->AddEntry(fits[i], names[i].c_str(), "l");
        //fit_legend->AddEntry(fits[i], Form("mean: %.2f/sigma: %.2f", fits[i]->GetParameter(1), fits[i]->GetParameter(2)), "");
    }
    fit_legend->Draw();

    TLatex latex_seed;
    latex_seed.SetTextSize(0.04);
    latex_seed.SetTextFont(62);
    latex_seed.DrawLatexNDC(0.63, 0.88, "Seed charge");
    latex_seed.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("p%s/%s/Thd%se", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));
    if(model == "masetti") {
        latex_seed.DrawLatexNDC(0.63, 0.73, "Masetti-Canali");
    }
    if(model == "jacoboni") {
        latex_seed.DrawLatexNDC(0.63, 0.73, "Jacoboni-Canali");
    }
    if(model == "canali") {
        latex_seed.DrawLatexNDC(0.63, 0.73, "Canali");
    }

    output_ROOTFile->cd();
    output_ROOTFile->cd(model.c_str());
    c_seed->Write(Form("hSeedCharge_%s_%s_Thd%s", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));

    //Fitting function
    //vFitSeedCharge = plot_MobilityModel::get_landauFit(vSeedCharge);
}

void plot_ComparingVoltage::get_voltageResidual(std::vector<TFile*> input_ROOTFile, TFile* output_ROOTFile, std::vector<std::string> voltage) {
    const std::string hResidual_dir_path = "DetectorHistogrammer/CE65/residuals/";

    const std::vector<std::string> axis = {"x", "y"};

    std::vector<TH1D*> vResidualX = {};
    std::vector<TH1D*> vResidualY = {};
    std::vector<TF1*> vFitResidualX = {};
    std::vector<TF1*> vFitResidualY = {};

    for(int i=0; i<input_ROOTFile.size(); ++i) {
        TH1D* hResidualX = (TH1D*)input_ROOTFile[i]->Get("DetectorHistogrammer/CE65/residuals/residual_x");
        if(hResidualX == nullptr) {
            LOG_WARNING.source("plot_ComparingVoltage::get_voltageResidual") << "Histogram is nullptr.";
            break;
        }
        vResidualX.push_back(hResidualX);
    }

    for(int i=0; i<input_ROOTFile.size(); ++i) {
        TH1D* hResidualY = (TH1D*)input_ROOTFile[i]->Get("DetectorHistogrammer/CE65/residuals/residual_y");
        if(hResidualY == nullptr) {
            LOG_WARNING.source("plot_ComparingVoltage::get_voltageResidual") << "Histogram is nullptr.";
            break;
        }
        vResidualY.push_back(hResidualY);
    }
    
    TCanvas* c_x = new TCanvas("c_x", "c_x", 800, 600);
    c_x->SetTopMargin(0.062);
    c_x->SetBottomMargin(0.14);
    c_x->SetLeftMargin(0.10);
    c_x->SetRightMargin(0.03);

    LOG_DEBUG.source("plot_ComparingVoltage::get_voltageResidual") << "Start drawing hist";
    for(int i=0; i<vResidualX.size(); ++i) {
        if(!vResidualX[i]) {
            LOG_WARNING.source("plot_ComparingVoltage::get_voltageResidual") << "vResidualX[" << i << "] is nullptr! Skipping this entry.";
            continue;
        }
        vResidualX[i]->SetMarkerSize(1.5);
        vResidualX[i]->SetMarkerStyle(20);
        if(i == 0) {
            vResidualX[i]->SetTitle(";x_{MC} - x_{cluster}; events");
            vResidualX[i]->GetXaxis()->SetRangeUser(x_Residual_min, x_Residual_max);
            vResidualX[i]->GetYaxis()->SetRangeUser(y_Residual_min, y_Residual_max);
            vResidualX[i]->GetXaxis()->SetTitleOffset(0.8);
            vResidualX[i]->GetYaxis()->SetTitleOffset(0.8);
            vResidualX[i]->GetXaxis()->SetTitleSize(0.06);
            vResidualX[i]->GetYaxis()->SetTitleSize(0.06);
            vResidualX[i]->GetXaxis()->SetLabelSize(0.04);
            vResidualX[i]->GetYaxis()->SetLabelSize(0.04);
            vResidualX[i]->SetMarkerColor(1);
            vResidualX[i]->Draw("PE");
        } else {
            vResidualX[i]->SetMarkerColor(i+1);
            vResidualX[i]->SetLineColor(i+1);
            if(i == 2) {
                vResidualX[i]->SetLineColor(kGreen+2);
                vResidualX[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                vResidualX[i]->SetLineColor(kCyan);
                vResidualX[i]->SetMarkerColor(kCyan);
            }
            vResidualX[i]->Draw("samePE");
        }
    }

    vFitResidualX = plot_MobilityModel::get_fits(vResidualX);

    // Write TLegend on TCanvas
    const double y_legend_max = 0.72 - (0.06*vResidualX.size());
    const double y_fit_legend_max = 0.88 - (0.071*vResidualX.size());
    TLegend* legend = new TLegend(0.7, 0.72, 0.88, y_legend_max);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    for(int i=0; i<vResidualX.size(); ++i) {
        legend->AddEntry(vResidualX[i], Form("%sV", voltage[i].c_str()), "p");
    }
    legend->Draw();
    TLegend* fit_legend = new TLegend(0.12, 0.88, 0.4, y_fit_legend_max);
    fit_legend->SetFillStyle(0);
    fit_legend->SetTextSize(0.03);
    fit_legend->SetBorderSize(0);
    for(int i=0; i<vResidualX.size(); ++i) {
        fit_legend->AddEntry(vFitResidualX[i], Form("%sV", voltage[i].c_str()), "");
        fit_legend->AddEntry(vFitResidualX[i], Form("mean: %.2f/sigma: %.2f", vFitResidualX[i]->GetParameter(1), vFitResidualX[i]->GetParameter(2)), "");
    }
    fit_legend->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatexNDC(0.63, 0.88, "Residual:x_{MC} - x_{cluster}");
    latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("p%s/%s/Thd%se", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));
    if(model == "masetti") {
        latex.DrawLatexNDC(0.63, 0.73, "Masetti-Canali");
    }
    if(model == "jacoboni") {
        latex.DrawLatexNDC(0.63, 0.73, "Jacoboni-Canali");
    }
    if(model == "canali") {
        latex.DrawLatexNDC(0.63, 0.73, "Canali");
    }

    output_ROOTFile->cd();
    output_ROOTFile->cd(model.c_str());
    c_x->Write(Form("hResidual_x_%s_%s_Thd%s", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));

    //
    TCanvas* c_y = new TCanvas("c_y", "c_y", 800, 600);
    c_y->SetTopMargin(0.062);
    c_y->SetBottomMargin(0.14);
    c_y->SetLeftMargin(0.10);
    c_y->SetRightMargin(0.03);

    for(int i=0; i<vResidualX.size(); ++i) {
        vResidualY[i]->SetMarkerSize(1.5);
        vResidualY[i]->SetMarkerStyle(20);
        if(i == 0) {
            vResidualY[i]->SetTitle(";y_{MC} - y_{cluster}; events");
            vResidualY[i]->GetXaxis()->SetRangeUser(x_Residual_min, x_Residual_max);
            vResidualY[i]->GetYaxis()->SetRangeUser(y_Residual_min, y_Residual_max);
            vResidualY[i]->GetXaxis()->SetTitleOffset(0.8);
            vResidualY[i]->GetYaxis()->SetTitleOffset(0.8);
            vResidualY[i]->GetXaxis()->SetTitleSize(0.06);
            vResidualY[i]->GetYaxis()->SetTitleSize(0.06);
            vResidualY[i]->GetXaxis()->SetLabelSize(0.04);
            vResidualY[i]->GetYaxis()->SetLabelSize(0.04);
            vResidualY[i]->SetMarkerColor(1);
            vResidualY[i]->Draw("PE");
        } else {
            vResidualY[i]->SetMarkerColor(i+1);
            vResidualY[i]->SetLineColor(i+1);
            if(i == 2) {
                vResidualY[i]->SetLineColor(kGreen+2);
                vResidualY[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                vResidualY[i]->SetLineColor(kCyan);
                vResidualY[i]->SetMarkerColor(kCyan);
            }
            vResidualY[i]->Draw("samePE");
        }
    }
    vFitResidualY = plot_MobilityModel::get_fits(vResidualY);

    // Write TLegend on TCanvas
    const double y_legend_y_max = 0.72 - (0.06*vResidualY.size());
    const double y_fit_legend_y_max = 0.88 - (0.071*vResidualY.size());
    TLegend* legend_y = new TLegend(0.7, 0.72, 0.88, y_legend_y_max);
    legend_y->SetFillStyle(0);
    legend_y->SetTextSize(0.04);
    legend_y->SetBorderSize(0);
    for(int i=0; i<vResidualY.size(); ++i) {
        legend_y->AddEntry(vResidualY[i], Form("%sV", voltage[i].c_str()), "p");
    }
    legend_y->Draw();
    TLegend* fit_legend_y = new TLegend(0.12, 0.88, 0.4, y_fit_legend_y_max);
    fit_legend_y->SetFillStyle(0);
    fit_legend_y->SetTextSize(0.03);
    fit_legend_y->SetBorderSize(0);
    for(int i=0; i<vResidualY.size(); ++i) {
        fit_legend_y->AddEntry(vFitResidualY[i], Form("%sV", voltage[i].c_str()), "");
        fit_legend_y->AddEntry(vFitResidualY[i], Form("mean: %.2f/sigma: %.2f", vFitResidualY[i]->GetParameter(1), vFitResidualY[i]->GetParameter(2)), "");
    }
    fit_legend_y->Draw();

    TLatex latex_y;
    latex_y.SetTextSize(0.04);
    latex_y.SetTextFont(62);
    latex_y.DrawLatexNDC(0.63, 0.88, "Residual:y_{MC} - y_{cluster}");
    latex_y.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("p%s/%s/Thd%se", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));
    if(model == "masetti") {
        latex_y.DrawLatexNDC(0.63, 0.73, "Masetti-Canali");
    }
    if(model == "jacoboni") {
        latex_y.DrawLatexNDC(0.63, 0.73, "Jacoboni-Canali");
    }
    if(model == "canali") {
        latex_y.DrawLatexNDC(0.63, 0.73, "Canali");
    }

    output_ROOTFile->cd();
    output_ROOTFile->cd(model.c_str());
    c_y->Write(Form("hResidual_y_%s_%s_Thd%s", pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str()));

}

void plot_ComparingVoltage::voltage_run() {
    set_rootStyle();

    const std::string hist_master_dir = "/home/towa/alice3/hist/ce65_sim_202505/";
    std::vector<TFile*> input_ROOTFile = {};

    std::string output_file_name = Form("plot/ce65_voltage.root");
    TFile* output = TFile::Open(output_file_name.c_str(), "RECREATE");
    TDirectory* masseti = output->mkdir(model_[0].c_str());
    TDirectory* canali = output->mkdir(model_[1].c_str());
    TDirectory* jacoboni = output->mkdir(model_[2].c_str());

    //model = "masetti";
    // pixel_pitch = pixel_pitch_[0];
    // chip_type = chip_type_[0];
    // threshold = threshold_[0];

    // plot_ComparingVoltage::get_voltageClsize(input_ROOTFile, output, voltage_);
    // plot_ComparingVoltage::get_voltageCharge(input_ROOTFile, output, voltage_);
    // plot_ComparingVoltage::get_voltageResidual(input_ROOTFile, output, voltage_);

    for(std::string model_loop_var : model_) {
        ::model = model_loop_var;
        for(std::string pixel_pitch_loop_var : pixel_pitch_) {
            ::pixel_pitch = pixel_pitch_loop_var;
            for(std::string chip_type_loop_var : chip_type_) {
                ::chip_type = chip_type_loop_var;
                for(std::string threshold_loop_var : threshold_) {
                    ::threshold = threshold_loop_var;
                    input_ROOTFile = {};
                    for(std::string voltage_loop_var : voltage_) {
                        ::voltage = voltage_loop_var;
                        TFile* input_file = TFile::Open(Form("%sn%sv/ce65_p%s_%s_Thd%se_%s.root", hist_master_dir.c_str(), voltage.c_str(), pixel_pitch.c_str(), chip_type.c_str(), threshold.c_str(), model.c_str()));
                        input_ROOTFile.push_back(input_file);
                    }

                    plot_ComparingVoltage::get_voltageClsize(input_ROOTFile, output, voltage_);
                    plot_ComparingVoltage::get_voltageCharge(input_ROOTFile, output, voltage_);
                    plot_ComparingVoltage::get_voltageResidual(input_ROOTFile, output, voltage_);
                }
            }
        }
    }

    output->Close();

    // Convert to PDF
    for(std::string model_loop_var : model_) {
        ::model = model_loop_var;
        plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), model.c_str(), Form("plot/ce65_voltage_1chip_%s.pdf", model.c_str()));
    } 
}
