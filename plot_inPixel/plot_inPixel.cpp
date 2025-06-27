#include "plot_inPixel.h"

plot_inPixel::plot_inPixel() {
    PIXEL_PITCH_    = {"15", "22p5"};
    CHIP_TYPE_      = {"std", "gap"};
    VOLTAGE_        = {"4", "7", "10"};
    MODEL_          = {"masetti"};
    THRESHOLD_      = {"0", "10", "25", "50", "100", "200", "300", "400"};
    //PIXEL_NUMBER_ = {"one_pixel", "multi_pixel"};
    PIXEL_NUMBER_ = {"multi_pixel"};
    TIME_ = plot_histogram::currentDateTime();
}

void plot_inPixel::set_2DSURFStyle(TCanvas* canvas, TH2D* hist) {
    canvas->Clear();
    canvas->SetTopMargin(0.062);
    canvas->SetBottomMargin(0.1);
    canvas->SetLeftMargin(0.14);
    canvas->SetRightMargin(0.04);

    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetZaxis()->SetTitleOffset(0.8);

    hist->Draw("SURF1");\
    // gPad->Update();
    // palette = (TPaletteAxis*)gPad->GetPrimitive("palette");
    // if(!palette) {
    //     LOG_WARNING.source("plot_inPixel::set_2DSURFStyle") << "Cannot find TPaletteAxis object.";
    // }
    // palette->SetX1NDC(0.87);
    // palette->SetY1NDC(0.2);
    // palette->SetX2NDC(0.90);
    // palette->SetY2NDC(0.85);
    // palette->SetTitleOffset(0.8);
    // hist->GetZaxis()->SetTitle("");
    // TLatex paletteTitle;
    // //paletteTitle->SetNDC();
    // paletteTitle.SetTextFont(hist->GetZaxis()->GetTitleFont());
    // paletteTitle.SetTextSize(hist->GetZaxis()->GetTitleSize());
    // paletteTitle.SetTextAngle(90);
    // paletteTitle.SetTextAlign(22);
    
    // double posX = palette->GetX2NDC() + 0.05;
    // double posY = (palette->GetY1NDC() + palette->GetY2NDC()) / 2.0;

    // paletteTitle.DrawLatexNDC(posX, posY, colorTitle.c_str());
}

void plot_inPixel::run_inPixel() {
    LOG_STATUS.source("plot_inPixel::run_inPixel") << "Start run for in-pixel analysis.";
    gStyle->SetPalette(55);

    std::string output_file_name = "/home/towa/alice3/plotter/plot/inPixel.root";
    std::string driftFile_dir_path = "/home/towa/alice3/hist/ce65_sim_driftTime_202506/";
    TFile* output = TFile::Open(output_file_name.c_str(), "RECREATE");

    TDirectory* masetti = output->mkdir("masetti");
    TDirectory* jacoboni = output->mkdir("jacoboni");
    TDirectory* canali = output->mkdir("canali");

    std::vector<TDirectory*> outputDir = {masetti, canali, jacoboni};

    TFile* inputROOTFile = nullptr;
    TProfile2D* p_clSize = nullptr;
    TProfile2D* p_residual = nullptr;
    TProfile2D* p_driftTime = nullptr;
    TPaletteAxis* palette = nullptr;
    TH2D* hClSize = nullptr;
    TH2D* hResidual = nullptr;
    TH2D* hDriftTime = nullptr;
    TCanvas* canvas = new TCanvas("canvas","canvas",800, 600);
    TLatex title;
    title.SetTextSize(0.04);
    title.SetTextFont(62);
    TLatex condition;
    condition.SetTextSize(0.03);
    condition.SetTextFont(62);
    // canvas->SetTopMargin(0.062);
    // canvas->SetBottomMargin(0.14);
    // canvas->SetLeftMargin(0.11);
    // canvas->SetRightMargin(0.03);
    std::string chip_variation;
    std::string chip_variation_text;
    for(int i=0; i<PIXEL_NUMBER_.size(); i++) {
        for(int j=0; j<MODEL_.size(); j++) {
            for(int k=0; k<CHIP_TYPE_.size(); k++) {
                for(int l=0; l<PIXEL_PITCH_.size(); l++) {
                    for(int n=0; n<VOLTAGE_.size(); n++) {
                        for(int m=0; m<THRESHOLD_.size(); m++) {
                            inputROOTFile = TFile::Open(Form("%sn%sv/%s/ce65_p%s_%s_%sV_Thd%se_%s_0.root", driftFile_dir_path.c_str(), VOLTAGE_[n].c_str(), PIXEL_NUMBER_[i].c_str(), PIXEL_PITCH_[l].c_str(), CHIP_TYPE_[k].c_str(), VOLTAGE_[n].c_str(), THRESHOLD_[m].c_str(), MODEL_[j].c_str()));
                            chip_variation = Form("%s_%s_%sV_Thd%se_%s_%s", PIXEL_PITCH_[l].c_str(), CHIP_TYPE_[k].c_str(), VOLTAGE_[n].c_str(), THRESHOLD_[m].c_str(), MODEL_[j].c_str(), PIXEL_NUMBER_[i].c_str());
                            //chip_variation_text = replace_underscore_to_slash(chip_variation);
                            chip_variation_text = Form("p%s/%s/%sV/Thd%se/%s/%s", PIXEL_PITCH_[l].c_str(), CHIP_TYPE_[k].c_str(), VOLTAGE_[n].c_str(), THRESHOLD_[m].c_str(), MODEL_[j].c_str(), PIXEL_NUMBER_[i].c_str());
                            gStyle->SetPalette(kViridis);
                            // gStyle->SetTextFont(132);
                            // gStyle->SetTitleFont(132, "XYZ");
                            // gStyle->SetLabelFont(132, "XYZ");

                            p_clSize = (TProfile2D*)inputROOTFile->Get("AnalysisPixel/CE65/inPixel_cluster_size");
                            p_residual = (TProfile2D*)inputROOTFile->Get("AnalysisPixel/CE65/inPixel_residual");
                            p_driftTime = (TProfile2D*)inputROOTFile->Get("AnalysisPixel/CE65/inPixel_electron_driftTime");
                            
                            hClSize = convert_toTH2D(p_clSize);
                            hResidual = convert_toTH2D(p_residual);
                            hDriftTime = convert_toTH2D(p_driftTime);

                            hClSize->SetMinimum(0);
                            hClSize->SetMaximum(8);
                            hClSize->SetTitle(";x w/in pixel [um];y w/in pixel [um];cluster size");
                            plot_inPixel::set_2DSURFStyle(canvas, hClSize);
                            title.DrawLatexNDC(0.20, 0.89, "Mean of Cluster Size");
                            // hClSize->SetContour(100);
                            // hClSize->GetXaxis()->SetRangeUser(-15,15);
                            // hClSize->GetYaxis()->SetRangeUser(-15, 15);

                            condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                            condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                            condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                            outputDir[j]->cd();
                            canvas->Write(Form("clSize_ce65_%s", chip_variation.c_str()));
                            canvas->SaveAs(Form("inPixel/clSize_ce65_%s.pdf", chip_variation.c_str()));

                            hResidual->SetMinimum(0);
                            hResidual->SetMaximum(15);
                            hResidual->SetTitle(";x w/in pixel [um];y w/in pixel [um];r_{MC}-r_{hit} [um]");
                            set_2DSURFStyle(canvas, hResidual);   
                            title.DrawLatexNDC(0.20, 0.89, "Mean of Residual (r_{MC} - r_{hit})");
                            condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                            condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                            condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                            outputDir[j]->cd();
                            canvas->Write(Form("residual_ce65_%s", chip_variation.c_str()));
                            canvas->SaveAs(Form("inPixel/residual_ce65_%s.pdf", chip_variation.c_str()));

                            hDriftTime->SetMinimum(0);
                            hDriftTime->SetMaximum(20);
                            hDriftTime->SetTitle(";x w/in pixel [um];y w/in pixel [um];drift time [ns]");
                            set_2DSURFStyle(canvas, hDriftTime);
                            title.DrawLatexNDC(0.20, 0.89, "90% Charge Collection Time (electron)");
                            condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                            condition.DrawLatexNDC(0.20, 0.82, Form("Plotted on %s", TIME_.c_str()));
                            condition.DrawLatexNDC(0.20, 0.79, chip_variation_text.c_str());
                            outputDir[j]->cd();
                            canvas->Write(Form("driftTime_ce65_%s", chip_variation.c_str()));
                            canvas->SaveAs(Form("inPixel/driftTime_ce65_%s.pdf", chip_variation.c_str()));

                            inputROOTFile->Close();
                        } // THRESHOLD_
                    } // VOLTAGE_
                } // PIXEL_PITCH_
            } // CHIP_TYPE_
        } // MODEL_
    } // PIXEL_NUMBER_

    output->Close();

    for(int i=0; i<MODEL_.size(); i++) {
            plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), MODEL_[i].c_str(), Form("plot/ce65_inPixel_%s.pdf", MODEL_[i].c_str()));
        }
}

TH2D* plot_inPixel::convert_toTH2D(TProfile2D* profile2D) {
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

std::string plot_inPixel::replace_underscore_to_slash(std::string text) {
    std::replace(text.begin(), text.end(), '_', '/');
    return text;
}