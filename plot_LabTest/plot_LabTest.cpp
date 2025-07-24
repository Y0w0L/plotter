#include "plot_LabTest.h"

plot_LabTest::plot_LabTest() {
    FILENAME_STD_1V     = "/home/towa/its3/labtest/labtest_202408/data/STD_1V_Fe55/STD_1V_Fe55_merged_seedthr1000_nbhrthr300_1_3.root";
    FILENAME_STD_4V     = "/home/towa/its3/labtest/labtest_202408/data/STD_4V_Fe55/STD_4V_Fe55_merged_seedthr1000_nbhrthr300_1_3.root";
    FILENAME_STD_7V     = "/home/towa/its3/labtest/labtest_202408/data/STD_7V_Fe55/STD_7V_Fe55_merged_seedthr1000_nbhrthr300_1_4.root";
    FILENAME_STD_10V    = "/home/towa/its3/labtest/labtest_202408/data/STD_10V_Fe55/STD_10V_Fe55_merged_seedthr1000_nbhrthr300_1_4.root";
    FILENAME_BLK_1V     = "/home/towa/its3/labtest/labtest_202408/data/BLT_1V_Fe55/BLT_1V_Fe55_merged_seedthr1000_nbhrthr300_1_2.root";
    FILENAME_BLK_4V     = "/home/towa/its3/labtest/labtest_202408/data/BLT_4V_Fe55/BLT_4V_Fe55_merged_seedthr1000_nbhrthr300_1_2.root";
    FILENAME_BLK_7V     = "/home/towa/its3/labtest/labtest_202408/data/BLT_7V_Fe55/BLT_7V_Fe55_merged_seedthr1000_nbhrthr300_1_2.root";
    FILENAME_BLK_10V    = "/home/towa/its3/labtest/labtest_202408/data/BLT_10V_Fe55/BLT_10V_Fe55_merged_seedthr1000_nbhrthr300_1_8.root";
    FILENAME_GAP_1V     = "/home/towa/its3/labtest/labtest_202408/data/GAP_1V_Fe55/GAP_1V_Fe55_merged_seedthr1000_nbhrthr300_1_4.root";
    FILENAME_GAP_4V     = "/home/towa/its3/labtest/labtest_202408/data/GAP_4V_Fe55/GAP_4V_Fe55_merged_seedthr1000_nbhrthr300_11_14.root";
    FILENAME_GAP_7V     = "/home/towa/its3/labtest/labtest_202408/data/GAP_7V_Fe55/GAP_7V_Fe55_merged_seedthr1000_nbhrthr300_1_5.root";
    FILENAME_GAP_10V    = "/home/towa/its3/labtest/labtest_202408/data/GAP_10V_Fe55/GAP_10V_Fe55_merged_seedthr1000_nbhrthr300_1_7.root";
    FILENAMES_STD_ = {FILENAME_STD_1V, FILENAME_STD_4V, FILENAME_STD_7V, FILENAME_STD_10V};
    FILENAMES_BLK_ = {FILENAME_BLK_1V, FILENAME_BLK_4V, FILENAME_BLK_7V, FILENAME_BLK_10V};
    FILENAMES_GAP_ = {FILENAME_GAP_1V, FILENAME_GAP_4V, FILENAME_GAP_7V, FILENAME_GAP_10V};
    FILENAMES_ = {FILENAMES_STD_, FILENAMES_BLK_, FILENAMES_GAP_};
    
    CHIP_TYPE_ = {"STD", "BLK", "GAP"};
    VOLTAGE_ = {"1V", "4V", "7V", "10V"};

    TIME_ = plot_histogram::currentDateTime();
}

void plot_LabTest::set_TH1DStyle(TH1D* hist) {
    hist->SetMarkerSize(0.5);
    hist->GetXaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
}

void plot_LabTest::set_TH1DPosition(TCanvas* canvas) {
    canvas->Clear();
    canvas->SetTopMargin(0.062);
    canvas->SetBottomMargin(0.14);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.03);
}

void plot_LabTest::run_LabTest() {
    LOG_STATUS.source("plot_LabTest::run_LabTest") << "Start run for Lab test analysis.";

    TFile* output = TFile::Open("/home/towa/alice3/plotter/plot/labtest.root", "RECREATE");

    TDirectory* chip_type = output->mkdir("chip_type");
    TDirectory* voltage = output->mkdir("voltage");
    TDirectory* mean = output->mkdir("mean");

    TFile* inputROOTFile = nullptr;
    TH1D* h_CType = nullptr;
    TH1D* h_Voltage = nullptr;
    TF1* f_CType = nullptr;
    TF1* f_Voltage = nullptr;
    TCanvas* canvas = new TCanvas("canvas","canvas",800,600);
    TLatex title;
    title.SetTextSize(0.04);
    title.SetTextFont(62);
    TLatex condition;
    condition.SetTextSize(0.03);
    condition.SetTextFont(62);

    // chip type
    for(int i=0; i<VOLTAGE_.size(); ++i) {
        set_TH1DPosition(canvas);
        TLegend* legend = new TLegend(0.7, 0.76, 0.88, 0.88);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);
        legend->SetBorderSize(0);

        for(int j=0; j<CHIP_TYPE_.size(); ++j) {
            inputROOTFile = TFile::Open(FILENAMES_[j][i].c_str());
            h_CType = (TH1D*)inputROOTFile->Get("h_mxAmpAC_spectra[0]");
            LOG_STATUS.source("plot_LabTest::run_LabTest()") << VOLTAGE_[i] << CHIP_TYPE_[j];
            LOG_STATUS.source("plot_LabTest::run_LabTest()") << "Hist Name is " << h_CType->GetName() << ".";
            set_TH1DStyle(h_CType);
            h_CType->Rebin(16);
            h_CType->SetMarkerStyle(8);

            // // Fitting
            // if(j == 0) {
            //     TF1* fit = new TF1("f_0", "gaus", )
            // }

            if(j == 0) {
                TF1* fit = new TF1("f_0","gaus", 1000, 10000);
                fit->SetLineColor(1);
                h_CType->Fit(fit, "RLQ", "", 5000, 6500);
                h_CType->SetStats(0);
                h_CType->SetTitle(";signal [ADC];entries");
                h_CType->GetXaxis()->SetRangeUser(2000,9000);
                h_CType->GetYaxis()->SetRangeUser(0,5000);
                h_CType->SetMarkerColor(1);
                h_CType->SetLineColor(1);
                //h_CType->SetMarkerStyle(8);
                h_CType->Draw("PE");
            } else {
                h_CType->SetLineColor(j+1);
                h_CType->SetMarkerColor(j+1);
                if(j == 1) {
                    //h_CType->SetMarkerStyle(21);
                    TF1* fit = new TF1("f_1","gaus", 1000, 10000);
                    fit->SetLineColor(2);
                    h_CType->Fit(fit, "RLQ", "", 4000, 5000);
                }
                if(j == 2) {
                    h_CType->SetLineColor(kGreen+2);
                    h_CType->SetMarkerColor(kGreen+2);
                    //h_CType->SetMarkerStyle(22);
                    TF1* fit = new TF1("f_2","gaus", 1000, 10000);
                    fit->SetLineColor(2);
                    h_CType->Fit(fit, "RLQ", "", 4200, 5100);
                }
                if(j == 4) {
                    h_CType->SetLineColor(kCyan);
                    h_CType->SetMarkerColor(kCyan);
                }
                h_CType->Draw("samePE");
            }
            legend->AddEntry(h_CType, Form("%s/%s", CHIP_TYPE_[j].c_str(), VOLTAGE_[i].c_str()), "p");
        } // CHIP_TYPE_
        legend->Draw();

        chip_type->cd();
        canvas->Write(VOLTAGE_[i].c_str());
    } // VOLTAGE_

    output->Close();
}