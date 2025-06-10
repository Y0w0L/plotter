#include "plot_MobilityModel.h"

const int judge = 0;
const bool scaling_residual = false;
const bool scaling_clSize = true;
std::string voltage_model;

plot_MobilityModel::plot_MobilityModel()
    : input_fileDir(""), output_fileDir(""), output_file(""),
      pixel_pitchs({}), chip_types({}), chip_variations({}), threshold_values({}), thresholds({}), model_names({}), models({}),
      h_clSize(nullptr),
      clSize_x_max(0), clSize_x_max_highThd(0), charge_x_max(0),
      clSize_y_max(0), clSize_y_max_highThd(0), charge_y_max(0),
      x_residual_max(0), y_residual_max(0),
      y_inPixelResidual_max(0)
{
    LOG_INFO.source("plot_MobilityModel::plot_MobilityModel") << "Initialize plot_MobilityModel.";
}

void plot_MobilityModel::clSizeStyle(TH1D* hist, const int color, const int marker_style) {
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(1.5);
    //hist->SetLineWidth(1.0);
}

TH1D* plot_MobilityModel::get_hist(TFile* file, const std::string& hist_name) {
    if(file == nullptr) {
        LOG_WARNING.source("plot_MobilityModel::get_hist") << "TFile pointer is nullptr. Retrun value will be nullptr.";
        return nullptr;
    }

    // get TObject from file
    TObject* obj = file->Get(hist_name.c_str());
    if(obj == nullptr) {
        LOG_WARNING.source("plot_MobilityModel::get_hist") << "Object of " << hist_name << " is nullptr in " << file->GetName() << ". Return value will be nullptr.";
        return nullptr;
    }

    TH1D* hist = dynamic_cast<TH1D*>(obj);
    // found object but it isn't TH1D type
    if(hist == nullptr) {
        LOG_WARNING.source("plot_MobilityModel::get_hist") << "TH1D object of " << hist_name << " is nullptr in " << file->GetName() << ". Return value will be nullptr.";
        return nullptr;
    }

    return hist;
}

std::vector<TFile*> plot_MobilityModel::get_TFiles(const std::vector<std::string>& filenames) {
    std::vector<TFile*> files;
    files.reserve(filenames.size());

    for(const std::string& filename : filenames) {
        TFile* file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie()) {
            LOG_ERROR.source("plot_MobilityModel::get_TFiles") << "Cannot open file: " << filename << ". Return value will be nullptr";
            delete file;
            files.push_back(nullptr);
        } else {
            files.push_back(file);
        }
    }
    return files;
}

std::vector<TH1D*> plot_MobilityModel::get_hists(
    std::vector<TFile*>& files,
    const std::string& hist_name) {
    
    std::vector<TH1D*> hists;
    hists.reserve(files.size());

    for(TFile*& file : files) {
        TH1D* hist = plot_MobilityModel::get_hist(file, hist_name);
        hists.push_back(hist);
    }

    return hists;
}

std::vector<TF1*> plot_MobilityModel::get_fits(std::vector<TH1D*> hists) {
    LOG_DEBUG.source("plot_MobilityModel::get_fits") << "Start get_fits.";
    std::vector<TF1*> fits;
    fits.reserve(hists.size());

    for(int i=0; i<hists.size(); ++i) {
        if(hists[i] == nullptr) {
            LOG_WARNING.source("plot_MobilityModel::get_fits") << "Histogram is nullptr. Return value will be nullptr.";
            fits.push_back(nullptr);
            continue;
        }

        if(scaling_residual & plot_histogram::search_word(hists[i]->GetName(), "residual")) {
            double max_amplitude = hists[i]->GetMaximum();
            hists[i]->Scale(1/max_amplitude);
        }

        if(i == 0) {
            TF1* fit = plot_histogram::optimise_hist_gaus(hists[i], 1);
            fits.push_back(fit);
        }
        if(i == 2) {
            TF1* fit = plot_histogram::optimise_hist_gaus(hists[i], kGreen+2);
            fits.push_back(fit);
        }
        if(i == 4) {
            TF1* fit = plot_histogram::optimise_hist_gaus(hists[i], kCyan);
            fits.push_back(fit);
        } else {
            TF1* fit = plot_histogram::optimise_hist_gaus(hists[i], i+1);
            fits.push_back(fit);
        }
    }

    return fits;
}

std::vector<TF1*> plot_MobilityModel::get_landauFits(std::vector<TH1D*> hists) {
    LOG_DEBUG.source("plot_MobilityModel::get_landauFits") << "Start get_landauFits.";
    std::vector<TF1*> fits;
    fits.reserve(hists.size());

    for(int i=0; i<hists.size(); ++i) {
        if(hists[i] == nullptr) {
            LOG_WARNING.source("plot_MobilityModel::get_landauFits") << "Histogram is nullptr. Return value will be nullptr.";
            fits.push_back(nullptr);
            continue;
        }

        if(i == 0) {
            TF1* fit = plot_histogram::optimise_hist_langaus(hists[i], 1);
            fits.push_back(fit);
        }
        if(i == 2) {
            TF1* fit = plot_histogram::optimise_hist_langaus(hists[i], kGreen+2);
            fits.push_back(fit);
        }
        if(i == 4) {
            TF1* fit = plot_histogram::optimise_hist_langaus(hists[i], kCyan);
            fits.push_back(fit);
        } else {
            TF1* fit = plot_histogram::optimise_hist_langaus(hists[i], i+1);
            fits.push_back(fit);
        }
    }
    return fits;
}

void plot_MobilityModel::write_clSizePlots(std::vector<TH1D*> hists, const std::vector<std::string>& names, const int& x_max, const int& y_max) {
    LOG_DEBUG.source("plot_MobilityModel::write_clSizePlots") << "Write cluster size plots.";

    for(int i=0; i<hists.size(); ++i) {
        if(hists[i] == nullptr) {
            continue;
        }
        LOG_INFO.source("plot_MobilityModel.write_clSizePlots") << "Hist name is " << hists[i]->GetName() <<".";
        if(scaling_clSize) {
            hists[i]->Scale(1/hists[i]->Integral());
        }
        if(i == 0) {
            clSizeStyle(hists[i], 1, 20);
            hists[i]->SetTitle(";cluster size; clusters");
            if(scaling_clSize) {
                hists[i]->SetTitle(";cluster size; clusters (normalized)");
            }
            hists[i]->GetXaxis()->SetRangeUser(1, x_max);
            hists[i]->GetYaxis()->SetRangeUser(0, y_max);
            hists[i]->GetXaxis()->SetTitleOffset(0.8);
            hists[i]->GetYaxis()->SetTitleOffset(0.8);
            hists[i]->GetXaxis()->SetTitleSize(0.06);
            hists[i]->GetYaxis()->SetTitleSize(0.06);
            hists[i]->GetXaxis()->SetLabelSize(0.04);
            hists[i]->GetYaxis()->SetLabelSize(0.04);
            hists[i]->Draw("histE");
        } else {
            clSizeStyle(hists[i], i+1, 20);
            if(i == 2) {
                hists[i]->SetLineColor(kGreen+2);
                hists[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                hists[i]->SetLineColor(kCyan);
                hists[i]->SetMarkerColor(kCyan);
            }
            hists[i]->Draw("samehistE");
        }
        double y_legend_size = 0.76 - (0.05*names.size());
        TLegend* legend = new TLegend(0.7, 0.76, 0.88, y_legend_size);
        //TLegend* legend = new TLegend(0.60, 0.65, 0.88, 0.88);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->SetBorderSize(0);
        for(int i=0; i<names.size(); ++i) {
            legend->AddEntry(hists[i], names[i].c_str(), "p");
        }
        legend->Draw();
    }
}

std::vector<TH1D*> plot_MobilityModel::write_chargePlots(std::vector<TH1D*> hists, std::vector<std::string> names, int x_max, int y_max) {
    LOG_DEBUG.source("plot_MobilityModel::write_chargePlots") << "Write charge plots.";
    //std::vector<TH1D*> scaledHists;
    for(int i=0; i<hists.size(); ++i) {
        int nBins = hists[i]->GetNbinsX();
        double xMax = hists[i]->GetXaxis()->GetBinCenter(nBins)*1000;
        TH1D* cloneHist = (TH1D*)hists[i]->Clone();
        hists[i] = new TH1D(hists[i]->GetName(), hists[i]->GetTitle(), nBins, 0, xMax);
        for(int j=1; j<=nBins; ++j) {
            hists[i]->SetBinContent(j, cloneHist->GetBinContent(j));
        } 
        LOG_INFO.source("plot_MobilityModel::write_chargePlots") << "Hist name is " << hists[i]->GetName() << ".";
        if(i == 0) {
            clSizeStyle(hists[i], 1, 20);
            hists[i]->SetTitle(";charge [e]; events");
            hists[i]->GetXaxis()->SetRangeUser(0, x_max);
            if(judge == 1) {
                hists[i]->GetYaxis()->SetRangeUser(1, y_max);
            }
            else {
                hists[i]->GetYaxis()->SetRangeUser(0, y_max);
            }
            hists[i]->GetXaxis()->SetTitleOffset(0.8);
            hists[i]->GetYaxis()->SetTitleOffset(0.8);
            hists[i]->GetXaxis()->SetTitleSize(0.06);
            hists[i]->GetYaxis()->SetTitleSize(0.06);
            hists[i]->GetXaxis()->SetLabelSize(0.04);
            hists[i]->GetYaxis()->SetLabelSize(0.04);
            hists[i]->Draw("PE");
        } else {
            clSizeStyle(hists[i], i+1, 20);
            if(i == 2) {
                hists[i]->SetLineColor(kGreen+2);
                hists[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                hists[i]->SetLineColor(kCyan);
                hists[i]->SetMarkerColor(kCyan);
            }
            hists[i]->Draw("samePE");
        }
        double y_legend_size = 0.76 - (0.05*names.size());\
        double y_fit_legend_size = 0.88 - (0.071*names.size());
        TLegend* legend = new TLegend(0.7, 0.76, 0.88, y_legend_size);
        //TLegend* legend = new TLegend(0.60, 0.65, 0.88, 0.88);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->SetBorderSize(0);
        for(int i=0; i<names.size(); ++i) {
            legend->AddEntry(hists[i], names[i].c_str(), "p");
        }
        legend->Draw();
        TLegend* fit_legend = new TLegend(0.12, 0.88, 0.4, y_fit_legend_size);
        fit_legend->SetFillStyle(0);
        fit_legend->SetTextSize(0.03);
        fit_legend->SetBorderSize(0);
        for(int i=0; i<names.size(); ++i) {
            //fit_legend->AddEntry(fits[i], names[i].c_str(), "l");
            //fit_legend->AddEntry(fits[i], Form("mean: %.2f/sigma: %.2f", fits[i]->GetParameter(1), fits[i]->GetParameter(2)), "");
        }
        fit_legend->Draw();
    }
    return hists;
}

void plot_MobilityModel::write_residualPlots(std::vector<TH1D*> hists, std::vector<TF1*> fits, std::vector<std::string> names, int x_max, int y_max, std::string axis) {
    LOG_DEBUG.source("plot_MobilityModel::write_residualPlots") << "Write residual plots.";
    for(int i=0; i<hists.size(); ++i) {
        // if(scaling_residual) {
        //     double max_amplitude = hists[i]->GetMaximum();
        //     hists[i]->Scale(1/max_amplitude);
        //     plot_histogram::optimise_hist_gaus(hists[i], i+1);
        //     //fits[i]->Scale(1/max_amplitude);
        // }
        LOG_INFO.source("plot_MobilityModel::write_residualPlots") << "Hist name is " << hists[i]->GetName() << ".";
        if(i == 0) {
            clSizeStyle(hists[i], 1, 20);
            hists[i]->SetTitle(Form(";%s_{MC} - %s_{cluster} [um]; events", axis.c_str(), axis.c_str()));
            hists[i]->GetXaxis()->SetRangeUser(-x_max, x_max);
            hists[i]->GetYaxis()->SetRangeUser(0, y_max);
            hists[i]->GetXaxis()->SetTitleOffset(0.8);
            hists[i]->GetYaxis()->SetTitleOffset(0.8);
            hists[i]->GetXaxis()->SetTitleSize(0.06);
            hists[i]->GetYaxis()->SetTitleSize(0.06);
            hists[i]->GetXaxis()->SetLabelSize(0.04);
            hists[i]->GetYaxis()->SetLabelSize(0.04);
            hists[i]->Draw("PE");
        } else {
            clSizeStyle(hists[i], i+1, 20);
            if(i == 2) {
                hists[i]->SetLineColor(kGreen+2);
                hists[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                hists[i]->SetLineColor(kCyan);
                hists[i]->SetMarkerColor(kCyan);
            }
            hists[i]->Draw("samePE");
        }
        double y_legend_size = 0.76 - (0.05*names.size());
        double y_fit_legend_size = 0.88 - (0.071*names.size());
        TLegend* legend = new TLegend(0.7, 0.76, 0.88, y_legend_size);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->SetBorderSize(0);
        for(int i=0; i<names.size(); ++i) {
            legend->AddEntry(hists[i], names[i].c_str(), "p");
        }
        legend->Draw();
        TLegend* fit_legend = new TLegend(0.12, 0.88, 0.4, y_fit_legend_size);
        fit_legend->SetFillStyle(0);
        fit_legend->SetTextSize(0.03);
        fit_legend->SetBorderSize(0);
        for(int i=0; i<names.size(); ++i) {
            fit_legend->AddEntry(fits[i], names[i].c_str(), "");
            fit_legend->AddEntry(fits[i], Form("mean: %.2f/sigma: %.2f", fits[i]->GetParameter(1), fits[i]->GetParameter(2)), "");
        }
        fit_legend->Draw();
    }
}

void plot_MobilityModel::write_inPixelResidualPlots(std::vector<TH1D*> hists, std::vector<std::string> names, int x_max, int y_max, std::string axis) {
    LOG_DEBUG.source("plot_MobilityModel::write_inPixelResidualPlots") << "Write in-pixel residual plots.";
    double y_legend_size = 0;
    for(int i=0; i<hists.size(); ++i) {
        LOG_INFO.source("plot_MobilityModel::write_inPixelResidualPlots") << "Hist name is " << hists[i]->GetName() << ".";
        if(i==0) {
            clSizeStyle(hists[i], 1, 20);
            if(axis == "xx") {
                hists[i]->SetTitle(";x in-pixel[um]; MAD(x_{MC}-x_{cluster})[um]");
            }
            if(axis == "xy") {
                hists[i]->SetTitle(";x in-pixel[um]; MAD(y_{MC}-y_{cluster})[um]");
            }
            if(axis == "yx") {
                hists[i]->SetTitle(";y in-pixel[um]; MAD(x_{MC}-x_{cluster})[um]");
            }
            if(axis == "yy") {
                hists[i]->SetTitle(";y in-pixel[um]; MAD(y_{MC}-y_{cluster})[um]");
            }
            //hists[i]->GetXaxis()->SetRangeUser(, x_max);
            hists[i]->GetYaxis()->SetRangeUser(0, y_max);
            hists[i]->GetXaxis()->SetTitleOffset(0.8);
            hists[i]->GetYaxis()->SetTitleOffset(0.65);
            hists[i]->GetXaxis()->SetTitleSize(0.06);
            hists[i]->GetYaxis()->SetTitleSize(0.06);
            hists[i]->GetXaxis()->SetLabelSize(0.04);
            hists[i]->GetYaxis()->SetLabelSize(0.04);
            hists[i]->SetTitleFont(42, "XYZ");
            hists[i]->Draw("PE");
        } else {
            clSizeStyle(hists[i], i+1, 20);
            if(i == 2) {
                hists[i]->SetLineColor(kGreen+2);
                hists[i]->SetMarkerColor(kGreen+2);
            }
            if(i == 4) {
                hists[i]->SetLineColor(kCyan);
                hists[i]->SetMarkerColor(kCyan);
            }
            hists[i]->Draw("samePE");
        }
        y_legend_size = 0.76 - (0.05*names.size());
        TLegend* legend = new TLegend(0.7, 0.76, 0.88, y_legend_size);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->SetBorderSize(0);
        for(int i=0; i<names.size(); ++i) {
            legend->AddEntry(hists[i], names[i].c_str(), "p");
        }
        legend->Draw();
    }
}

void plot_MobilityModel::save_clSizePlots(TFile* output_file, std::vector<std::string> filenames, std::vector<std::string> names, std::string save_name, int x_max, int y_max) {
    std::vector<TFile*> files = get_TFiles(filenames);
    std::vector<TH1D*> hists = get_hists(files, "DetectorHistogrammer/CE65/cluster_size/cluster_size");

    for(int i=0; i<hists.size(); ++i) {
        if(hists[i] == nullptr) {
            LOG_WARNING.source("plot_MobilityModel::save_clSizePlots") << "Histogram is nullptr";
        }
    }

    // for(int i=0; i<filenames.size(); i++) {
    //     std::cout << filenames[i] << std::endl;
    // }
    // for(int i=0; i<names.size(); i++) {
    //     std::cout << names[i] << std::endl;
    // }

    std::string output_dir = "cluster_size/";

    LOG_DEBUG.source("plot_MobilityModel::save_clSizePlots") << "Save cluster size plots.";
    TCanvas* canvas = new TCanvas("c", "c", 800, 600);
    canvas->SetTopMargin(0.062);
    canvas->SetBottomMargin(0.14);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.03);
    write_clSizePlots(hists, names, x_max, y_max);
    // for(int i=0; i<files.size(); i++) {
    //     files[i]->Close();
    // }S
    std::string name = save_name;
    plot_histogram::replace_string(name, "_", "/");
    // TLegend* legend = new TLegend(0.55, 0.4, 0.88, 0.6);
    // legend->SetFillStyle(0);
    // legend->SetTextSize(0.04);
    // legend->SetBorderSize();
    // legend->SetHeader(name.c_str());
    // legend->Draw();
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatexNDC(0.63, 0.88, "Cluster size");
    latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("%s/%sV",name.c_str(),voltage_model.c_str()));

    output_file->cd();
    output_file->cd(output_dir.c_str());
    canvas->Write(save_name.c_str());
    // canvas->SaveAs("test.pdf");
    delete canvas;
}

void plot_MobilityModel::save_chargePlots(TFile* output_file, std::vector<std::string> filenames, std::vector<std::string> names, std::string save_name, int x_max, int y_max) {
    std::vector<TFile*> files = get_TFiles(filenames);
    std::vector<TH1D*> hists = get_hists(files, "DetectorHistogrammer/CE65/charge/cluster_charge");
    std::vector<TH1D*> scaledHists;
    scaledHists.reserve(hists.size());

    for(int i=0; i<hists.size(); ++i) {
        if(hists[i] == nullptr) {
            LOG_WARNING.source("plot_MobilityModel::save_chargePlots") << "Histogram is nullptr";
        }
    }

    // std::vector<TF1*> fits = get_landauFits(hists);

    // for(int i=0; i<filenames.size(); i++) {
    //     std::cout << filenames[i] << std::endl;
    // }
    // for(int i=0; i<names.size(); i++) {
    //     std::cout << names[i] << std::endl;
    // }

    std::string output_dir = "cluster_charge/";

    LOG_DEBUG.source("plot_MobilityModel::save_chargePlots") << "Save charge plots.";
    TCanvas* canvas = new TCanvas("c", "c", 800, 600);
    if(judge == 1) {
        canvas->SetLogy();
    }
    canvas->SetTopMargin(0.065);
    canvas->SetBottomMargin(0.12);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.03);
    scaledHists = write_chargePlots(hists, names, x_max, y_max);
    // for(int i=0; i<files.size(); i++) {
    //     files[i]->Close();
    // }

    std::vector<TF1*> fits = get_landauFits(scaledHists);

    std::string name = save_name;
    plot_histogram::replace_string(name, "_", "/");
    // TLegend* legend = new TLegend(0.55, 0.4, 0.88, 0.6);
    // legend->SetFillStyle(0);
    // legend->SetTextSize(0.04);
    // legend->SetBorderSize();
    // legend->SetHeader(name.c_str());
    // legend->Draw();
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatexNDC(0.63, 0.88, "Cluster charge");
    latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("%s/%sV",name.c_str(),voltage_model.c_str()));

    output_file->cd();
    output_file->cd(output_dir.c_str());
    canvas->Write(save_name.c_str());
    delete canvas;
}

void plot_MobilityModel::save_seedChargePlots(TFile* output_file, std::vector<std::string> filenames, std::vector<std::string> names, std::string save_name, int x_max, int y_max) {
    std::vector<TFile*> files = get_TFiles(filenames);
    std::vector<TH1D*> hists = get_hists(files, "DetectorHistogrammer/CE65/charge/seed_charge");
    std::vector<TH1D*> scaledHists;

    for(int i=0; i<hists.size(); ++i) {
        if(hists[i] == nullptr) {
            LOG_WARNING.source("plot_MobilityModel::save_seedChargePlots") << "Histogram is nullptr";
        }
    }

    // std::vector<TF1*> fits = get_landauFits(hists);

    // for(int i=0; i<filenames.size(); i++) {
    //     std::cout << filenames[i] << std::endl;
    // }
    // for(int i=0; i<names.size(); i++) {
    //     std::cout << names[i] << std::endl;
    // }

    std::string output_dir = "seed_charge/";

    TCanvas* canvas = new TCanvas("c", "c", 800, 600);
    if(judge == 1) {
        canvas->SetLogy();
    }
    canvas->SetTopMargin(0.065);
    canvas->SetBottomMargin(0.12);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.03);
    scaledHists = write_chargePlots(hists, names, x_max, y_max);
    // for(int i=0; i<files.size(); i++) {
    //     files[i]->Close();
    // }

    std::vector<TF1*> fits = get_landauFits(scaledHists);

    std::string name = save_name;
    plot_histogram::replace_string(name, "_", "/");
    // TLegend* legend = new TLegend(0.55, 0.4, 0.88, 0.6);
    // legend->SetFillStyle(0);
    // legend->SetTextSize(0.04);
    // legend->SetBorderSize();
    // legend->SetHeader(name.c_str());
    // legend->Draw();
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatexNDC(0.63, 0.88, "Seed charge");
    latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
    latex.DrawLatexNDC(0.63, 0.78, Form("%s/%sV",name.c_str(),voltage_model.c_str()));

    output_file->cd();
    output_file->cd(output_dir.c_str());
    canvas->Write(save_name.c_str());
    delete canvas;
}

void plot_MobilityModel::save_residualPlots(TFile* output_file, std::vector<std::string> filenames, std::vector<std::string> names, std::string save_name, int x_max, int y_max) {
    std::vector<TFile*> files = get_TFiles(filenames);
    std::vector<std::string> axes = {"x", "y"};
    for(auto axis : axes) {
        std::vector<TH1D*> hists = get_hists(files, Form("DetectorHistogrammer/CE65/residuals/residual_%s", axis.c_str()));
        std::vector<TF1*> fits = get_fits(hists);
        for(int i=0; i<hists.size(); ++i) {
            if(hists[i] == nullptr) {
                LOG_WARNING.source("plot_MobilityModel::save_residualPlots") << "Histogram is nullptr";
            }
        }
        std::string output_dir = "residuals/";
        TCanvas* canvas = new TCanvas("c", "c", 800, 600);
        canvas->SetTopMargin(0.065);
        canvas->SetBottomMargin(0.12);
        canvas->SetLeftMargin(0.10);
        canvas->SetRightMargin(0.03);

        write_residualPlots(hists, fits, names, x_max, y_max, axis);
    
        std::string name = save_name;
        plot_histogram::replace_string(name, "_", "/");
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetTextFont(62);
        latex.DrawLatexNDC(0.63, 0.88, Form("Residual:%s_{MC} - %s_{cluster}", axis.c_str(), axis.c_str()));
        latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
        latex.DrawLatexNDC(0.63, 0.78, Form("%s/%sV",name.c_str(),voltage_model.c_str()));

        // TLatex latex_files;
        // latex_files.SetTextSize(0.04);
        // latex_files.SetTextFont(62);
        // for(int i=0; i<filenames.size(); i++) {
        //     //std::cout << filenames[i] << std::endl;
        //     latex_files.DrawLatexNDC(0.8,1-0.05*i, filenames[i].c_str());
        // }

        output_file->cd();
        output_file->cd(output_dir.c_str());
        canvas->Write(save_name.c_str());
        delete canvas;
    }
}

void plot_MobilityModel::save_inPixelResidualPlots(TFile* output_file, std::vector<std::string> filenames, std::vector<std::string> names, std::string save_name, int x_max, int y_max) {
    std::vector<TFile*> files = get_TFiles(filenames);
    std::vector<std::string> axes = {"xx", "xy", "yx", "yy"}; // "xx" x_vs_x (1st_vs_2nd)
    std::vector<TH1D*> hists;
    std::string output_dir = "inPixelResiduals/";
    for(int i=0; i<axes.size(); ++i) {
        if(axes[i] == "xx") {
            hists = get_hists(files, "DetectorHistogrammer/CE65/residuals/residual_x_vs_x");
        }
        if(axes[i] == "xy") {
            hists = get_hists(files, "DetectorHistogrammer/CE65/residuals/residual_x_vs_y");
        }
        if(axes[i] == "yx") {
            hists = get_hists(files, "DetectorHistogrammer/CE65/residuals/residual_y_vs_x");
        }
        if(axes[i] == "yy") {
            hists = get_hists(files, "DetectorHistogrammer/CE65/residuals/residual_y_vs_y");
        }
        for(int i=0; i<hists.size(); ++i) {
            if(hists[i] == nullptr) {
                LOG_WARNING.source("plot_MobilityModel::save_inPixelResidualPlots") << "Histogram is nullptr.";
            }
        }
        TCanvas* canvas = new TCanvas("c", "c", 800, 600);
        canvas->SetTopMargin(0.065);
        canvas->SetBottomMargin(0.12);
        canvas->SetLeftMargin(0.10);
        canvas->SetRightMargin(0.03);

        write_inPixelResidualPlots(hists, names, x_max, y_max, axes[i]);

        std::string name = save_name;
        plot_histogram::replace_string(name, "_", "/");
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetTextFont(62);
        latex.DrawLatexNDC(0.63, 0.88, "In-pixel residual");
        latex.DrawLatexNDC(0.63, 0.83, "Electron:3GeV/c");
        latex.DrawLatexNDC(0.63, 0.78, Form("%s/%sV",name.c_str(),voltage_model.c_str()));

        output_file->cd();
        output_file->cd(output_dir.c_str());
        canvas->Write(save_name.c_str());
        delete canvas;
    }
}
// TODO need to be fixed 
void plot_MobilityModel::save_suminPixelResidualPlots(TFile* output_file, std::vector<std::string> filenames, std::vector<std::string> names, std::string save_name, int x_max, int y_max) {
    LOG_WARNING.source("plot_MobilityModel::save_suminPixelResidualPlots") << "Save sum in-pixel residual plots";
    std::vector<TFile*> files = get_TFiles(filenames);
    std::vector<std::string> axes = {"xx", "xy", "yx", "yy"}; // "xx" x_vs_x (1st_vs_2nd)
    std::vector<TH1D*> hists;
    std::string output_dir = "inPixelSumResiduals/";
    TCanvas* canvas = new TCanvas("c", "c", 800, 600);
    canvas->SetTopMargin(0.065);
    canvas->SetBottomMargin(0.12);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.03);

    for(int i=0; i<files.size(); ++i) {
        std::vector<TH1D*> inPixel_residuals_x = {};
        std::vector<TH1D*> inPixel_residuals_y = {};
        LOG_INFO.source("plot_MobilityModel::save_suminPixelResidualPlots") << "File name is " << files[i]->GetName() << ".";
        for(int j=0; j<axes.size(); ++j) {
            if(axes[j] == "xx") {
                LOG_INFO.source("plot_MobilityModel::save_suminPixelResidualPlots") << "Get hist from DetectorHistogrammer/CE65/residuals/residual_x_vs_x";
                TH1D* hist_xx = get_hist(files[i], "DetectorHistogrammer/CE65/residuals/residual_x_vs_x");
                inPixel_residuals_x.push_back(hist_xx);
            }
            if(axes[j] == "xy") {
                LOG_INFO.source("plot_MobilityModel::save_suminPixelResidualPlots") << "Get hist from DetectorHistogrammer/CE65/residuals/residual_x_vs_y";
                TH1D* hist_xy = get_hist(files[i], "DetectorHistogrammer/CE65/residuals/residual_x_vs_y");
                inPixel_residuals_x.push_back(hist_xy);
            }
            // if(axes[j] == "yx") {
            //     std::cout << "Get hist: " << "DetectorHistogrammer/CE65/residuals/residual_y_vs_x" << std::endl;
            //     TH1D* hist_yx = get_hist(files[i], "DetectorHistogrammer/CE65/residuals/residual_y_vs_x");
            //     inPixel_residuals_y.push_back(hist_yx);
            // }
            // if(axes[j] == "yy") {
            //     std::cout << "Get hist: " << "DetectorHistogrammer/CE65/residuals/residual_y_vs_y" << std::endl;
            //     TH1D* hist_yy = get_hist(files[i], "DetectorHistogrammer/CE65/residuals/residual_y_vs_y");
            //     inPixel_residuals_y.push_back(hist_yy);
            // }
        }
        for(int j=0; j<inPixel_residuals_x.size(); ++j) {
            if(inPixel_residuals_x[j] == nullptr) {
                LOG_WARNING.source("plot_MobilityModel::save_suminPixelResidualPlots") << "Histogram is nullptr.";
                break;
            }
            LOG_INFO.source("plot_MobilityModel::save_suminPixelResidualPlots") << "Hist name -suminPixel: " << inPixel_residuals_x[j]->GetName() << ".";
        }
        // sum x
        TH1D* sum_residual_x = new TH1D("sum_residual_x", ";x in-pixel[um]; MAD(r)[um]", inPixel_residuals_x[0]->GetNbinsX(), inPixel_residuals_x[0]->GetXaxis()->GetXmin(), inPixel_residuals_x[0]->GetXaxis()->GetXmax());
        for(int j=1; j<= inPixel_residuals_x[0]->GetNbinsX(); ++j) {
            double y_xx = inPixel_residuals_x[0]->GetBinContent(j);
            double y_xy = inPixel_residuals_x[1]->GetBinContent(j);
            double y_sum = sqrt(y_xx*y_xx + y_xy*y_xy);
            sum_residual_x->SetBinContent(j, y_sum);
        }
        if(i == 0) {
            sum_residual_x->SetTitle(";x in-pixel[um]; MAD(r)[um]");
            sum_residual_x->SetLineColor(1);
            sum_residual_x->SetMarkerColor(1);
            sum_residual_x->GetXaxis()->SetTitleOffset(0.8);
            sum_residual_x->GetYaxis()->SetTitleOffset(0.65);
            sum_residual_x->GetXaxis()->SetTitleSize(0.06);
            sum_residual_x->GetYaxis()->SetTitleSize(0.06);
            sum_residual_x->GetXaxis()->SetLabelSize(0.04);
            sum_residual_x->GetYaxis()->SetLabelSize(0.04);
            sum_residual_x->SetTitleFont(42, "XYZ");
            sum_residual_x->GetYaxis()->SetRangeUser(0, 30);
            sum_residual_x->Draw("PE");
        } else {
            sum_residual_x->SetLineColor(i+1);
            sum_residual_x->SetMarkerColor(i+1);
            sum_residual_x->Draw("samePE");
        }
    }

    output_file->cd();
    output_file->cd(output_dir.c_str());
    canvas->Write(save_name.c_str());
    delete canvas;
}

// void plot_MobilityModel::all_plot_run() {
//     set_rootStyle();

//     std::cout << "Start run for all plots" << std::endl;
//     std::vector<std::string> chip_type = {"std", "gap"};
//     std::vector<std::string> pixel_pitch = {"p15", "p22p5"};
//     std::vector<std::string> voltage_model = {"4", "7", "10"};
//     std::vector<std::string> threshold = {"0", "10", "25", "50", "100", "200", "300", "400"};
//     //std::vector<std::string> model = {"canali", "jacoboni", "masetti_canali"};
//     std::vector<std::string> model = {"masetti_canali"};

//     std::string hist_master_dir_path = "/home/towa/alice3/hist/ce65_sim_202505/";

//     for (int i=0; i<voltage_model.size(); i++) {
//         std::string hist_path = form("%sn%sv/", hist_master_dir_path.c_str(), voltage_model[i].c_str());
//         for (int j=0; j<model.size(); j++) {
//             for (int k=0; k<chip_type.size(); k++) {
//                 for (int l=0; l<pixel_pitch.size(); l++) {
//                     for (int m=0; m<threshold.size(); m++) {
//                         std::string filename = form("%sce65_%s_%s_%s_%s_%s.root", hist_path.c_str(),pixel_pitch[l].c_str(), chip_type[k].c_str(), threshold[m].c_str(), model[j].c_str());
//                         std::cout << "File name: " filename << std::endl;
//                     }
//                 }
//             }
//         }
//     }
// }

void plot_MobilityModel::run() {
    LOG_STATUS.source("plot_MobilityModel::run") << "Starat Main process in plot_MobilityModel.";
    set_rootStyle();

    voltage_model = "4";

    if(judge == 0) {
        LOG_STATUS.source("plot_MobilityModel.run") << "Start run for 1chip simulation data";
        input_fileDir = Form("/home/towa/alice3/hist/ce65_sim_202505/n%sv/", voltage_model.c_str());
        output_file = Form("ce65_model_1chip_n%sv.root", voltage_model.c_str());
    }
    if(judge == 1) {
        LOG_STATUS.source("plot_MobilityModel.run") << "Start run for lab simulation data";
        input_fileDir = "/home/towa/alice3/hist/lab_sim/";
        output_file = "ce65_model_lab.root";
    }

    //input_fileDir = "/home/towa/alice3/hist/ce65_sim/";
    //output_fileDir = "/home/towa/alice3/hist/ce65_sim/";
    output_fileDir = "plot/";
    //output_file = "ce65_model.root";

    pixel_pitchs = {"15", "22p5"};
    chip_types = {"gap", "std"};
    for(const auto& chip : chip_types) {
        for(const auto& pitch : pixel_pitchs) {
            chip_variations.push_back("p" + pitch+ "_" + chip);
        }
    }
    threshold_values = {"0", "10", "25", "50", "100", "200"};
    for(const auto& th : threshold_values) {
        thresholds.push_back("Thd" + th + "e");
    }
    model_names = {"canali", "jacoboni", "masetti"};
    models = {"Canali", "Jacoboni-Canali", "Masetti-Canali"};

    std::sort(pixel_pitchs.begin(), pixel_pitchs.end());
    std::sort(chip_types.begin(), chip_types.end());
    std::sort(chip_variations.begin(), chip_variations.end());
    //std::sort(threshold_values.begin(), threshold_values.end());
    //std::sort(thresholds.begin(), thresholds.end());
    std::sort(model_names.begin(), model_names.end());
    std::sort(models.begin(), models.end());

    std::vector<std::string> filenames;

    output_file = output_fileDir + output_file;

    for (const auto& pitch : pixel_pitchs) {
        for (const auto& chip : chip_types) {
            for (const auto& th : thresholds) {
                for (const auto& model : model_names) {
                    std::string filename = input_fileDir + "ce65_p" + pitch + "_" + chip + "_" + th + "_" + model + ".root";
                    filenames.push_back(filename);
                }
            }
        }
    }

    LOG_STATUS.source("plot_MobilityModel::run") << "Start file filtering.";
    std::vector<std::string> filenames_p15_std_Thd0e = plot_histogram::and_filter_filenames(filenames, {"p15", "std","Thd0e"});
    std::vector<std::string> filenames_p15_gap_Thd0e = plot_histogram::and_filter_filenames(filenames, {"p15", "gap","Thd0e"});
    std::vector<std::string> filenames_p22p5_std_Thd0e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std","Thd0e"});
    std::vector<std::string> filenames_p22p5_gap_Thd0e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap","Thd0e"});
    std::vector<std::string> filenames_p15_std_Thd10e = plot_histogram::and_filter_filenames(filenames, {"p15", "std","Thd10e"});
    std::vector<std::string> filenames_p15_gap_Thd10e = plot_histogram::and_filter_filenames(filenames, {"p15", "gap","Thd10e"});
    std::vector<std::string> filenames_p22p5_std_Thd10e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std","Thd10e"});
    std::vector<std::string> filenames_p22p5_gap_Thd10e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap","Thd10e"});
    std::vector<std::string> filenames_p15_std_Thd25e = plot_histogram::and_filter_filenames(filenames, {"p15", "std","Thd25e"});
    std::vector<std::string> filenames_p15_gap_Thd25e = plot_histogram::and_filter_filenames(filenames, {"p15", "gap","Thd25e"});
    std::vector<std::string> filenames_p22p5_std_Thd25e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std","Thd25e"});
    std::vector<std::string> filenames_p22p5_gap_Thd25e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap","Thd25e"});
    std::vector<std::string> filenames_p15_std_Thd50e = plot_histogram::and_filter_filenames(filenames, {"p15", "std","Thd50e"});
    std::vector<std::string> filenames_p15_gap_Thd50e = plot_histogram::and_filter_filenames(filenames, {"p15", "gap","Thd50e"});
    std::vector<std::string> filenames_p22p5_std_Thd50e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std","Thd50e"});
    std::vector<std::string> filenames_p22p5_gap_Thd50e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap","Thd50e"});
    std::vector<std::string> filenames_p15_std_Thd100e = plot_histogram::and_filter_filenames(filenames, {"p15", "std","Thd100e"});
    std::vector<std::string> filenames_p15_gap_Thd100e = plot_histogram::and_filter_filenames(filenames, {"p15", "gap","Thd100e"});
    std::vector<std::string> filenames_p22p5_std_Thd100e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std","Thd100e"});
    std::vector<std::string> filenames_p22p5_gap_Thd100e = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap","Thd100e"});
    std::vector<std::string> filenames_Thd0e_canali = plot_histogram::and_filter_filenames(filenames, {"Thd0e", "canali."});
    std::vector<std::string> filenames_Thd0e_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"Thd0e", "jacoboni"});
    std::vector<std::string> filenames_Thd0e_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"Thd0e", "masetti"});
    std::vector<std::string> filenames_Thd10e_canali = plot_histogram::and_filter_filenames(filenames, {"Thd10e", "canali"});
    std::vector<std::string> filenames_Thd10e_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"Thd10e", "jacoboni"});
    std::vector<std::string> filenames_Thd10e_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"Thd10e", "masetti"});
    std::vector<std::string> filenames_Thd25e_canali = plot_histogram::and_filter_filenames(filenames, {"Thd25e", "canali"});
    std::vector<std::string> filenames_Thd25e_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"Thd25e", "jacoboni"});
    std::vector<std::string> filenames_Thd25e_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"Thd25e", "masetti"});
    std::vector<std::string> filenames_Thd50e_canali = plot_histogram::and_filter_filenames(filenames, {"Thd50e", "canali"});
    std::vector<std::string> filenames_Thd50e_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"Thd50e", "jacoboni"});
    std::vector<std::string> filenames_Thd50e_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"Thd50e", "masetti"});
    std::vector<std::string> filenames_Thd100e_canali = plot_histogram::and_filter_filenames(filenames, {"Thd100e", "canali"});
    std::vector<std::string> filenames_Thd100e_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"Thd100e", "jacoboni"});
    std::vector<std::string> filenames_Thd100e_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"Thd100e", "masetti"});
    std::vector<std::string> filenames_p15_std_canali = plot_histogram::and_filter_filenames(filenames, {"p15", "std", "canali"});
    std::vector<std::string> filenames_p15_gap_canali = plot_histogram::and_filter_filenames(filenames, {"p15", "gap", "canali"});
    std::vector<std::string> filenames_p22p5_std_canali = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std", "canali"});
    std::vector<std::string> filenames_p22p5_gap_canali = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap", "canali"});
    std::vector<std::string> filenames_p15_std_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"p15", "std", "jacoboni"});
    std::vector<std::string> filenames_p15_gap_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"p15", "gap", "jacoboni"});
    std::vector<std::string> filenames_p22p5_std_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std", "jacoboni"});
    std::vector<std::string> filenames_p22p5_gap_jacoboniCanali = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap", "jacoboni"});
    std::vector<std::string> filenames_p15_std_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"p15", "std", "masetti"});
    std::vector<std::string> filenames_p15_gap_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"p15", "gap", "masetti"});
    std::vector<std::string> filenames_p22p5_std_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"p22p5", "std", "masetti"});
    std::vector<std::string> filenames_p22p5_gap_masettiCanali = plot_histogram::and_filter_filenames(filenames, {"p22p5", "gap", "masetti"});
    
    TFile* output = TFile::Open(output_file.c_str(), "RECREATE");
    TDirectory* clSize = output->mkdir("cluster_size");
    TDirectory* cluster_charge = output->mkdir("cluster_charge");
    TDirectory* seed_charge = output->mkdir("seed_charge");
    TDirectory* residual = output->mkdir("residuals");
    TDirectory* inPixelResidual = output->mkdir("inPixelResiduals");
    TDirectory* inPixelSumResidual = output->mkdir("inPixelSumResiduals");

    // for 1chip & not scaled
    if(judge == 0) {
        clSize_x_max = 15;
        clSize_x_max_highThd = 8;
        charge_x_max = 3000;
        x_residual_max = 15;

        clSize_y_max = 800;
        charge_y_max = 700;
        y_residual_max = 800;
        y_inPixelResidual_max = 15;
        if(scaling_residual) {
            y_residual_max = 1.2;
        }
        if(scaling_clSize) {
            clSize_y_max = 0.8;
        }
    }
    // for lab test & not scaled
    if(judge == 1) {
        clSize_x_max = 15;
        clSize_x_max_highThd = 8;
        charge_x_max = 3000;
        x_residual_max = 20;

        clSize_y_max = 15000;
        charge_y_max = 20000;
        y_residual_max = 300;
    }

    // save cluster size plots
    save_clSizePlots(output, filenames_p15_std_Thd0e, models, "p15_std_Thd0e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_Thd0e, models, "p15_gap_Thd0e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_Thd0e, models, "p22p5_std_Thd0e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_Thd0e, models, "p22p5_gap_Thd0e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_Thd10e, models, "p15_std_Thd10e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_Thd10e, models, "p15_gap_Thd10e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_Thd10e, models, "p22p5_std_Thd10e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_Thd10e, models, "p22p5_gap_Thd10e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_Thd25e, models, "p15_std_Thd25e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_Thd25e, models, "p15_gap_Thd25e", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_Thd25e, models, "p22p5_std_Thd25e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_Thd25e, models, "p22p5_gap_Thd25e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_Thd50e, models, "p15_std_Thd50e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_Thd50e, models, "p15_gap_Thd50e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_Thd50e, models, "p22p5_std_Thd50e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_Thd50e, models, "p22p5_gap_Thd50e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_Thd100e, models, "p15_std_Thd100e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_Thd100e, models, "p15_gap_Thd100e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_Thd100e, models, "p22p5_std_Thd100e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_Thd100e, models, "p22p5_gap_Thd100e", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd0e_canali, chip_variations, "Thd0e_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd0e_jacoboniCanali, chip_variations, "Thd0e_Jacoboni-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd0e_masettiCanali, chip_variations, "Thd0e_Masetti-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd10e_canali, chip_variations, "Thd10e_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd10e_jacoboniCanali, chip_variations, "Thd10e_Jacoboni-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd10e_masettiCanali, chip_variations, "Thd10e_Masetti-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd25e_canali, chip_variations, "Thd25e_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd25e_jacoboniCanali, chip_variations, "Thd25e_Jacoboni-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd25e_masettiCanali, chip_variations, "Thd25e_Masetti-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd50e_canali, chip_variations, "Thd50e_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd50e_jacoboniCanali, chip_variations, "Thd50e_Jacoboni-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd50e_masettiCanali, chip_variations, "Thd50e_Masetti-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd100e_canali, chip_variations, "Thd100e_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd100e_jacoboniCanali, chip_variations, "Thd100e_Jacoboni-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_Thd100e_masettiCanali, chip_variations, "Thd100e_Masetti-Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_canali, thresholds, "p15_std_Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_canali, thresholds, "p15_gap_Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_canali, thresholds, "p22p5_std_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_canali, thresholds, "p22p5_gap_Canali", clSize_x_max_highThd, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_jacoboniCanali, thresholds, "p15_std_Jacoboni-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_jacoboniCanali, thresholds, "p15_gap_Jacoboni-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_jacoboniCanali, thresholds, "p22p5_std_Jacoboni-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_jacoboniCanali, thresholds, "p22p5_gap_Jacoboni-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_std_masettiCanali, thresholds, "p15_std_Masetti-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p15_gap_masettiCanali, thresholds, "p15_gap_Masetti-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_std_masettiCanali, thresholds, "p22p5_std_Masetti-Canali", clSize_x_max, clSize_y_max);
    save_clSizePlots(output, filenames_p22p5_gap_masettiCanali, thresholds, "p22p5_gap_Masetti-Canali", clSize_x_max, clSize_y_max);

    // save charge plots
    save_chargePlots(output, filenames_p15_std_Thd0e, models, "p15_std_Thd0e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_Thd0e, models, "p15_gap_Thd0e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_Thd0e, models, "p22p5_std_Thd0e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_Thd0e, models, "p22p5_gap_Thd0e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_Thd10e, models, "p15_std_Thd10e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_Thd10e, models, "p15_gap_Thd10e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_Thd10e, models, "p22p5_std_Thd10e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_Thd10e, models, "p22p5_gap_Thd10e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_Thd25e, models, "p15_std_Thd25e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_Thd25e, models, "p15_gap_Thd25e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_Thd25e, models, "p22p5_std_Thd25e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_Thd25e, models, "p22p5_gap_Thd25e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_Thd50e, models, "p15_std_Thd50e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_Thd50e, models, "p15_gap_Thd50e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_Thd50e, models, "p22p5_std_Thd50e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_Thd50e, models, "p22p5_gap_Thd50e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_Thd100e, models, "p15_std_Thd100e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_Thd100e, models, "p15_gap_Thd100e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_Thd100e, models, "p22p5_std_Thd100e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_Thd100e, models, "p22p5_gap_Thd100e", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd0e_canali, chip_variations, "Thd0e_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd0e_jacoboniCanali, chip_variations, "Thd0e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd0e_masettiCanali, chip_variations, "Thd0e_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd10e_canali, chip_variations, "Thd10e_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd10e_jacoboniCanali, chip_variations, "Thd10e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd10e_masettiCanali, chip_variations, "Thd10e_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd25e_canali, chip_variations, "Thd25e_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd25e_jacoboniCanali, chip_variations, "Thd25e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd25e_masettiCanali, chip_variations, "Thd25e_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd50e_canali, chip_variations, "Thd50e_Canali",charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd50e_jacoboniCanali, chip_variations, "Thd50e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd50e_masettiCanali, chip_variations, "Thd50e_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd100e_canali, chip_variations, "Thd100e_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd100e_jacoboniCanali, chip_variations, "Thd100e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_Thd100e_masettiCanali, chip_variations, "Thd100e_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_canali, thresholds, "p15_std_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_canali, thresholds, "p15_gap_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_canali, thresholds, "p22p5_std_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_canali, thresholds, "p22p5_gap_Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_jacoboniCanali, thresholds, "p15_std_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_jacoboniCanali, thresholds, "p15_gap_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_jacoboniCanali, thresholds, "p22p5_std_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_jacoboniCanali, thresholds, "p22p5_gap_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_std_masettiCanali, thresholds, "p15_std_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p15_gap_masettiCanali, thresholds, "p15_gap_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_std_masettiCanali, thresholds, "p22p5_std_Masetti-Canali", charge_x_max, charge_y_max);
    save_chargePlots(output, filenames_p22p5_gap_masettiCanali, thresholds, "p22p5_gap_Masetti-Canali", charge_x_max, charge_y_max);

    // save seed_charge plots
    save_seedChargePlots(output, filenames_p15_std_Thd0e, models, "p15_std_Thd0e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_Thd0e, models, "p15_gap_Thd0e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_Thd0e, models, "p22p5_std_Thd0e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_Thd0e, models, "p22p5_gap_Thd0e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_Thd10e, models, "p15_std_Thd10e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_Thd10e, models, "p15_gap_Thd10e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_Thd10e, models, "p22p5_std_Thd10e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_Thd10e, models, "p22p5_gap_Thd10e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_Thd25e, models, "p15_std_Thd25e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_Thd25e, models, "p15_gap_Thd25e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_Thd25e, models, "p22p5_std_Thd25e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_Thd25e, models, "p22p5_gap_Thd25e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_Thd50e, models, "p15_std_Thd50e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_Thd50e, models, "p15_gap_Thd50e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_Thd50e, models, "p22p5_std_Thd50e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_Thd50e, models, "p22p5_gap_Thd50e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_Thd100e, models, "p15_std_Thd100e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_Thd100e, models, "p15_gap_Thd100e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_Thd100e, models, "p22p5_std_Thd100e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_Thd100e, models, "p22p5_gap_Thd100e", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd0e_canali, chip_variations, "Thd0e_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd0e_jacoboniCanali, chip_variations, "Thd0e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd0e_masettiCanali, chip_variations, "Thd0e_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd10e_canali, chip_variations, "Thd10e_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd10e_jacoboniCanali, chip_variations, "Thd10e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd10e_masettiCanali, chip_variations, "Thd10e_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd25e_canali, chip_variations, "Thd25e_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd25e_jacoboniCanali, chip_variations, "Thd25e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd25e_masettiCanali, chip_variations, "Thd25e_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd50e_canali, chip_variations, "Thd50e_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd50e_jacoboniCanali, chip_variations, "Thd50e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd50e_masettiCanali, chip_variations, "Thd50e_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd100e_canali, chip_variations, "Thd100e_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd100e_jacoboniCanali, chip_variations, "Thd100e_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_Thd100e_masettiCanali, chip_variations, "Thd100e_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_canali, thresholds, "p15_std_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_canali, thresholds, "p15_gap_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_canali, thresholds, "p22p5_std_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_canali, thresholds, "p22p5_gap_Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_jacoboniCanali, thresholds, "p15_std_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_jacoboniCanali, thresholds, "p15_gap_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_jacoboniCanali, thresholds, "p22p5_std_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_jacoboniCanali, thresholds, "p22p5_gap_Jacoboni-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_std_masettiCanali, thresholds, "p15_std_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p15_gap_masettiCanali, thresholds, "p15_gap_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_std_masettiCanali, thresholds, "p22p5_std_Masetti-Canali", charge_x_max, charge_y_max);
    save_seedChargePlots(output, filenames_p22p5_gap_masettiCanali, thresholds, "p22p5_gap_Masetti-Canali", charge_x_max, charge_y_max);

    // save residual plots
    save_residualPlots(output, filenames_p15_std_Thd0e, models, "p15_std_Thd0e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_Thd0e, models, "p15_gap_Thd0e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_Thd0e, models, "p22p5_std_Thd0e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_Thd0e, models, "p22p5_gap_Thd0e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_Thd10e, models, "p15_std_Thd10e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_Thd10e, models, "p15_gap_Thd10e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_Thd10e, models, "p22p5_std_Thd10e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_Thd10e, models, "p22p5_gap_Thd10e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_Thd25e, models, "p15_std_Thd25e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_Thd25e, models, "p15_gap_Thd25e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_Thd25e, models, "p22p5_std_Thd25e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_Thd25e, models, "p22p5_gap_Thd25e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_Thd50e, models, "p15_std_Thd50e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_Thd50e, models, "p15_gap_Thd50e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_Thd50e, models, "p22p5_std_Thd50e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_Thd50e, models, "p22p5_gap_Thd50e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_Thd100e, models, "p15_std_Thd100e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_Thd100e, models, "p15_gap_Thd100e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_Thd100e, models, "p22p5_std_Thd100e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_Thd100e, models, "p22p5_gap_Thd100e", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd0e_canali, chip_variations, "Thd0e_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd0e_jacoboniCanali, chip_variations, "Thd0e_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd0e_masettiCanali, chip_variations, "Thd0e_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd10e_canali, chip_variations, "Thd10e_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd10e_jacoboniCanali, chip_variations, "Thd10e_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd10e_masettiCanali, chip_variations, "Thd10e_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd25e_canali, chip_variations, "Thd25e_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd25e_jacoboniCanali, chip_variations, "Thd25e_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd25e_masettiCanali, chip_variations, "Thd25e_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd50e_canali, chip_variations, "Thd50e_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd50e_jacoboniCanali, chip_variations, "Thd50e_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd50e_masettiCanali, chip_variations, "Thd50e_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd100e_canali, chip_variations, "Thd100e_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd100e_jacoboniCanali, chip_variations, "Thd100e_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_Thd100e_masettiCanali, chip_variations, "Thd100e_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_canali, thresholds, "p15_std_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_canali, thresholds, "p15_gap_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_canali, thresholds, "p22p5_std_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_canali, thresholds, "p22p5_gap_Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_jacoboniCanali, thresholds, "p15_std_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_jacoboniCanali, thresholds, "p15_gap_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_jacoboniCanali, thresholds, "p22p5_std_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_jacoboniCanali, thresholds, "p22p5_gap_Jacoboni-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_std_masettiCanali, thresholds, "p15_std_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p15_gap_masettiCanali, thresholds, "p15_gap_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_std_masettiCanali, thresholds, "p22p5_std_Masetti-Canali", x_residual_max, y_residual_max);
    save_residualPlots(output, filenames_p22p5_gap_masettiCanali, thresholds, "p22p5_gap_Masetti-Canali", x_residual_max, y_residual_max);

    save_inPixelResidualPlots(output, filenames_p15_std_Thd0e, models, "p15_std_Thd0e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_Thd0e, models, "p15_gap_Thd0e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_Thd0e, models, "p22p5_std_Thd0e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_Thd0e, models, "p22p5_gap_Thd0e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_Thd10e, models, "p15_std_Thd10e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_Thd10e, models, "p15_gap_Thd10e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_Thd10e, models, "p22p5_std_Thd10e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_Thd10e, models, "p22p5_gap_Thd10e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_Thd25e, models, "p15_std_Thd25e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_Thd25e, models, "p15_gap_Thd25e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_Thd25e, models, "p22p5_std_Thd25e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_Thd25e, models, "p22p5_gap_Thd25e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_Thd50e, models, "p15_std_Thd50e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_Thd50e, models, "p15_gap_Thd50e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_Thd50e, models, "p22p5_std_Thd50e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_Thd50e, models, "p22p5_gap_Thd50e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_Thd100e, models, "p15_std_Thd100e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_Thd100e, models, "p15_gap_Thd100e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_Thd100e, models, "p22p5_std_Thd100e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_Thd100e, models, "p22p5_gap_Thd100e", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd0e_canali, chip_variations, "Thd0e_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd0e_jacoboniCanali, chip_variations, "Thd0e_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd0e_masettiCanali, chip_variations, "Thd0e_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd10e_canali, chip_variations, "Thd10e_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd10e_jacoboniCanali, chip_variations, "Thd10e_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd10e_masettiCanali, chip_variations, "Thd10e_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd25e_canali, chip_variations, "Thd25e_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd25e_jacoboniCanali, chip_variations, "Thd25e_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd25e_masettiCanali, chip_variations, "Thd25e_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd50e_canali, chip_variations, "Thd50e_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd50e_jacoboniCanali, chip_variations, "Thd50e_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd50e_masettiCanali, chip_variations, "Thd50e_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd100e_canali, chip_variations, "Thd100e_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd100e_jacoboniCanali, chip_variations, "Thd100e_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_Thd100e_masettiCanali, chip_variations, "Thd100e_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_canali, thresholds, "p15_std_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_canali, thresholds, "p15_gap_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_canali, thresholds, "p22p5_std_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_canali, thresholds, "p22p5_gap_Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_jacoboniCanali, thresholds, "p15_std_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_jacoboniCanali, thresholds, "p15_gap_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_jacoboniCanali, thresholds, "p22p5_std_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_jacoboniCanali, thresholds, "p22p5_gap_Jacoboni-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_std_masettiCanali, thresholds, "p15_std_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p15_gap_masettiCanali, thresholds, "p15_gap_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_std_masettiCanali, thresholds, "p22p5_std_Masetti-Canali", x_residual_max, y_inPixelResidual_max);
    save_inPixelResidualPlots(output, filenames_p22p5_gap_masettiCanali, thresholds, "p22p5_gap_Masetti-Canali", x_residual_max, y_inPixelResidual_max);

    save_suminPixelResidualPlots(output, filenames_p15_std_Thd0e, models, "p15_std_Thd0e", x_residual_max, y_inPixelResidual_max);

    output->Close();
    LOG_STATUS.source("plot_MobilityModel::run") << "Finish writing histograms to " << output_file << ".";

    if(judge == 0) {
        plot_histogram::saveCanvasesToPDF(output_file.c_str(), "cluster_size", Form("plot/ce65_model_1chip_n%sv_clSize.pdf", voltage_model.c_str()));
        plot_histogram::saveCanvasesToPDF(output_file.c_str(), "cluster_charge", Form("plot/ce65_model_1chip_n%sv_clCharge.pdf", voltage_model.c_str()));
        plot_histogram::saveCanvasesToPDF(output_file.c_str(), "seed_charge", Form("plot/ce65_model_1chip_n%sv_seCharge.pdf", voltage_model.c_str()));
        plot_histogram::saveCanvasesToPDF(output_file.c_str(), "residuals", Form("plot/ce65_model_1chip_n%sv_residuals.pdf", voltage_model.c_str()));
        plot_histogram::saveCanvasesToPDF(output_file.c_str(), "inPixelResiduals", Form("plot/ce65_model_1chip_n%sv_inPixelResiduals.pdf", voltage_model.c_str()));
    }
    if(judge == 1) {
        plot_histogram::saveCanvasesToPDF("ce65_model_lab.root", "cluster_size", "ce65_model_lab_clSize.pdf");
        plot_histogram::saveCanvasesToPDF("ce65_model_lab.root", "charge", "ce65_model_lab_charge.pdf");
        plot_histogram::saveCanvasesToPDF("ce65_model_lab.root", "residuals", "ce65_model_lab_residuals.pdf");
    }

    // TCanvas* c = new TCanvas("c", "c", 800, 600);
    // TH1D* hist = get_hist(TFile::Open("/home/towa/alice3/hist/ce65_sim/ce65_p15_gap_Thd0e_jacoboniCanali.root"), "DetectorHistogrammer/CE65/residuals/residual_x");
    // TF1* doublegaus = plot_histogram::optimise_hist_doublegaus(hist, 1);
    // hist->Draw();
    // hist->Fit(doublegaus, "RL");
    // c->SaveAs("residual_x.pdf");

    
    // std::vector<std::string> model_of_propagation = {"masetti_canali", "canali", "jacoboni_canali"};
    // std::vector<std::string> chip_type = {"std", "gap"};
    // std::vector<std::string> pixel_pitch = {"15", "22p5"};
    // std::vector<std::string> voltage_model = {"4", "7", "10"};
    // std::vector<std::string> threshold = {"0", "10", "25", "50", "100", "200", "300", "400"};

    // hist_dir_path = "/home/towa/alice3/hist/ce65_sim_202505/"


}