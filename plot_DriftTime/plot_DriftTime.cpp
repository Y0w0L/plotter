#include "plot_DriftTime.h"

const std::vector<std::string> CHIP_TYPE_ = {"std", "gap"};
const std::vector<std::string> PIXEL_PITCH_ = {"15", "22p5"};
const std::vector<std::string> VOLTAGE_ = {"4", "7", "10"};
const std::vector<std::string> MODEL_ = {"masetti", "canali", "jacoboni"};
const std::vector<std::string> SOURCE_POSITION_NAME_ = {"0", "1", "2", "3", "4", "5"};
const std::vector<std::string> PIXEL_NUMBER_ = {"one_pixel", "multi_pixel"};

const std::string TIME_ = plot_histogram::currentDateTime();

plot_DriftTime::plot_DriftTime() {
    LOG_STATUS.source("plot_DriftTime::plot_DriftTime") << "plot_DriftTime object is created.";
}

TCanvas* plot_DriftTime::get_driftTime(std::vector<TFile*> inputROOTFile, TFile* outputROOTFile, const std::vector<std::string>& position_name) {
    LOG_DEBUG.source("plot_DriftTime::get_driftTime") << "Start get_driftTime.";
    const std::string hDriftTime_path = "GenericPropagation/CE65/drift_time_histo";
    std::vector<TH1D*> vDriftTime = {};

    const int size_ofVector = position_name.size();
    vDriftTime.reserve(size_ofVector);
    TH1D* hClSize = (TH1D*)inputROOTFile[0]->Get("DetectorHistogrammer/CE65/cluster_size/cluster_size");
    //double scaleFactor = hClSize->GetEntries();
    TH1D* hDriftTime = nullptr;
    for(int i=0; i<size_ofVector; i++) {
        hDriftTime = (TH1D*)inputROOTFile[i]->Get(hDriftTime_path.c_str());
        if(hDriftTime == nullptr) {
            LOG_WARNING.source("plot_DriftTime::get_driftTime") << "Histogram is nullptr.";
            break;
        }
        //hDriftTime->Scale(1.0/5.0);
        vDriftTime.push_back(hDriftTime);
    }
    
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    c->SetTopMargin(0.062);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.03);

    TLegend* legend = new TLegend(0.50, 0.50, 0.89, 0.92);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);

    for(int i=0; i<size_ofVector; i++) {
        if(!vDriftTime[i]) {
            LOG_WARNING.source("plot_DriftTime::get_driftTime") << "histogram is nullptr. Skip writing this histogram on canvas.";
            continue;
        }
        vDriftTime[i]->SetMarkerStyle(20);
        vDriftTime[i]->SetMarkerSize(1.0);
        vDriftTime[i]->SetLineWidth(1); 
        if(i == 0) {
            vDriftTime[i]->SetTitle(";drift time[ns];charge carriers");
            vDriftTime[i]->GetXaxis()->SetRangeUser(0, 25);
            vDriftTime[i]->GetYaxis()->SetRangeUser(0, 1000);
            vDriftTime[i]->GetXaxis()->SetTitleOffset(0.8);
            vDriftTime[i]->GetYaxis()->SetTitleOffset(0.8);
            vDriftTime[i]->GetXaxis()->SetTitleSize(0.06);
            vDriftTime[i]->GetYaxis()->SetTitleSize(0.06);
            vDriftTime[i]->GetXaxis()->SetLabelSize(0.04);
            vDriftTime[i]->GetYaxis()->SetLabelSize(0.04);
            vDriftTime[i]->SetMarkerColor(1);
            vDriftTime[i]->SetLineColor(1);
            vDriftTime[i]->Draw("LP");
        }
        else {
            if(i == 2) {
                vDriftTime[i]->SetMarkerColor(kGreen+2);
                vDriftTime[i]->SetLineColor(kGreen+2);
            }
            if(i == 4) {
                vDriftTime[i]->SetMarkerColor(kCyan);
                vDriftTime[i]->SetLineColor(kCyan);
            } else {
                vDriftTime[i]->SetMarkerColor(i+1);
                vDriftTime[i]->SetLineColor(i+1);
            }
            vDriftTime[i]->Draw("sameLP");
        }
    }
    for(int i=0; i<size_ofVector; i++) {
        legend->AddEntry(vDriftTime[i], Form("position: %s/ Mean: %.2f", SOURCE_POSITION_NAME_[i].c_str(), vDriftTime[i]->GetMean()), "p");
    }
    legend->Draw();
    //c->SaveAs("driftTime.root");
    // outputROOTFile->cd();
    // c->Write("driftTime");
    //c->SetLogy();
    return c;
}

void plot_DriftTime::run_driftTime() {
    LOG_STATUS.source("plot_DriftTime::run_driftTime") << "Start run for DriftTime.";
    std::vector<TFile*> inputROOTFile;
    std::string driftFile_dir_path = "/home/towa/alice3/hist/ce65_sim_driftTime_202506/";
    std::string output_file_name = "/home/towa/alice3/plotter/plot/driftTime.root";
    TFile* output = TFile::Open(output_file_name.c_str(), "RECREATE");

    TDirectory* masetti = output->mkdir("masetti");
    TDirectory* jacoboni = output->mkdir("jacoboni");
    TDirectory* canali = output->mkdir("canali");

    std::vector<TDirectory*> outputDir = {masetti, canali, jacoboni};
    
    for(int m=0; m<PIXEL_NUMBER_.size();m++) {
        for(int i=0; i<MODEL_.size(); i++) {
            for(int j=0; j<CHIP_TYPE_.size(); j++) {
                for(int k=0; k<PIXEL_PITCH_.size(); k++) {
                    for(int l=0; l<VOLTAGE_.size(); l++) {
                        inputROOTFile = {};
                        for(int n=0; n<SOURCE_POSITION_NAME_.size(); n++) {
                            TFile* input_file = TFile::Open(Form("%sn%sv/%s/ce65_p%s_%s_%sV_%s_%s.root", driftFile_dir_path.c_str(), VOLTAGE_[l].c_str(), PIXEL_NUMBER_[m].c_str(), PIXEL_PITCH_[k].c_str(), CHIP_TYPE_[j].c_str(), VOLTAGE_[l].c_str(), MODEL_[i].c_str(), SOURCE_POSITION_NAME_[n].c_str()));
                            inputROOTFile.push_back(input_file);
                        } // SOURCE_POSITION_NAME_
                        TCanvas* c = get_driftTime(inputROOTFile, output, SOURCE_POSITION_NAME_);
                        TLatex title;
                        title.SetTextSize(0.04);
                        title.SetTextFont(62);
                        title.DrawLatexNDC(0.20, 0.88, "Drift time");
                        TLatex condition;
                        condition.SetTextSize(0.02);
                        condition.SetTextFont(62);
                        condition.DrawLatexNDC(0.20, 0.85, "Electron:3GeV/c");
                        condition.DrawLatexNDC(0.20, 0.83, Form("Plotted on %s", TIME_.c_str()));
                        condition.DrawLatexNDC(0.20, 0.81, Form("%s/%s/p%s/%sV/%s", PIXEL_NUMBER_[m].c_str(), CHIP_TYPE_[j].c_str(), PIXEL_PITCH_[k].c_str(), VOLTAGE_[l].c_str(), MODEL_[i].c_str()));
                        output->cd();
                        outputDir[i]->cd();
                        c->Write(Form("driftTime_ce65_%s_%s_%sV_%s_%s", CHIP_TYPE_[j].c_str(), PIXEL_PITCH_[k].c_str(), VOLTAGE_[l].c_str(), MODEL_[i].c_str(), PIXEL_NUMBER_[m].c_str()));
                    } // VOLTAGE_
                } // PIXEL_PITCH_
            } // CHIP_TYPE_
        } // MODEL_
    } // PIXEL_NUMBER_
        
    LOG_STATUS.source("plot_DriftTime::run_driftTime") << "make histogram";

    // for(int i=0; i<MODEL_.size(); i++) {
    //     for(int j=0; j<CHIP_TYPE_.size(); j++) {
    //         for(int k=0; k<PIXEL_PITCH_.size(); k++) {
    //             for(int l=0; l<VOLTAGE_.size(); l++) {
    //                 get_driftTime(inputROOTFile, output, SOURCE_POSITION_NAME_);
    //             } // VOLTAGE_
    //         } // PIXEL_PITCH_
    //     } // CHIP_TYPE_
    // } // MODEL_

    output->Close();
    for(int i=0; i<MODEL_.size(); i++) {
            plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), MODEL_[i].c_str(), Form("plot/ce65_driftTime_%s.pdf", MODEL_[i].c_str()));
        }
    // for(int j=0; j<PIXEL_NUMBER_.size(); j++) {
    //     for(int i=0; i<MODEL_.size(); i++) {
    //         plot_histogram::saveCanvasesToPDF(output_file_name.c_str(), MODEL_[i].c_str(), Form("plot/ce65_driftTime_%s_%s.pdf", PIXEL_NUMBER_[j].c_str(), MODEL_[i].c_str()));
    //     }
    // }
}