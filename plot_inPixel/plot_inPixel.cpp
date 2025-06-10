#include "plot_inPixel.h"

plot_inPixel::plot_inPixel() {
    PIXEL_PITCH_    = {"15", "22p5"};
    CHIP_TYPE_      = {"std", "gap"};
    VOLTAGE_        = {"4", "7", "10"};
    MODEL_          = {"masetti"};
    THRESHOLD_      = {"0", "10", "25", "50", "100", "200", "300", "400"};
    //PIXEL_NUMBER_ = {"one_pixel", "multi_pixel"};
    PIXEL_NUMBER_ = {"multi_pixel"};
}

void plot_inPixel::run_inPixel() {
    LOG_STATUS.source("plot_inPixel::run_inPixel") << "Start run for in-pixel analysis.";
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
    TH2D* hClSize = nullptr;
    TH2D* hResidual = nullptr;
    TH2D* hDriftTime = nullptr;
    for(int i=0; i<PIXEL_NUMBER_.size(); i++) {
        for(int j=0; j<MODEL_.size(); j++) {
            for(int k=0; k<CHIP_TYPE_.size(); k++) {
                for(int l=0; l<PIXEL_PITCH_.size(); l++) {
                    for(int n=0; n<VOLTAGE_.size(); n++) {
                        for(int m=0; m<THRESHOLD_.size(); m++) {
                            inputROOTFile = TFile::Open(Form("%sn%sv/%s/ce65_p%s_%s_%sV_Thd%se_%s_0.root", driftFile_dir_path.c_str(), VOLTAGE_[n].c_str(), PIXEL_NUMBER_[i].c_str(), PIXEL_PITCH_[l].c_str(), CHIP_TYPE_[k].c_str(), VOLTAGE_[n].c_str(), THRESHOLD_[m].c_str(), MODEL_[j].c_str()));
                            p_clSize = (TProfile2D*)inputROOTFile->Get("AnalysisPixel/CE65/inPixel_cluster_size");
                            p_residual = (TProfile2D*)inputROOTFile->Get("AnalysisPixel/CE65/inPixel_residual");
                            p_driftTime = (TProfile2D*)inputROOTFile->Get("AnalysisPixel/CE65/inPixel_electron_dirftTime");
                            
                            hClSize = convert_toTH2D(p_clSize);
                            hResidual = convert_toTH2D(p_residual);
                            hDriftTime = convert_toTH2D(p_driftTime);

                            TCanvas* c_clSize = new TCanvas("c_clSize", "c_clSize", 800, 600);
                            c_clSize->SetTopMargin(0.062);
                            c_clSize->SetBottomMargin(0.14);
                            c_clSize->SetLeftMargin(0.11);
                            c_clSize->SetRightMargin(0.03);

                            hClSize->Draw("COLZ");
                            outputDir[j]->cd();
                            output->Write();

                        } // THRESHOLD_
                    } // VOLTAGE_
                } // PIXEL_PITCH_
            } // CHIP_TYPE_
        } // MODEL_
    } // PIXEL_NUMBER_

    output->Close();
}

TH2D* plot_inPixel::convert_toTH2D(TProfile2D* profile2D) {
    if(!profile2D) {
        return nullptr;
    }

    // get axis information from TProfile2D
    int nbinsX  = profile2D->GetNbinsX();
    double xMin = profile2D->GetXaxis()->GetXmin();
    double xMax = profile2D->GetXaxis()->GetXmax();

    int nbinsY  = profile2D->GetNbinsX();
    double yMin = profile2D->GetYaxis()->GetXmin();
    double yMax = profile2D->GetYaxis()->GetXmax();

    // make new TH2D
    TH2D* h2d = new TH2D(profile2D->GetName(), profile2D->GetTitle(), nbinsX, xMin, xMax, nbinsY, yMin, yMax);

    double content;
    double error;
    // copy contents and fill
    for(int i=1; i<=nbinsX; i++) {
        for(int j=1; j<nbinsY; j++) {
            content = profile2D->GetBinContent(i,j);
            error = profile2D->GetBinError(i,j);

            h2d->SetBinContent(i,j,content);
            h2d->SetBinError(i,j,error);
        }
    }

    return h2d;
}