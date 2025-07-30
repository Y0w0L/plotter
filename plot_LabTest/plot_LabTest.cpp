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
    
    CHIP_TYPE_ = {"std", "blk", "gap"};
    VOLTAGE_ = {"1", "4", "7", "10"};

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

    TCanvas* canvas = new TCanvas("canvas","canvas",800,600);
    TLatex title;
    TLatex condition;
    title.SetTextSize(0.04);
    title.SetTextFont(62);
    condition.SetTextSize(0.03);
    condition.SetTextFont(62);

    bool charge_as_chipType = true;
    if(charge_as_chipType) {
        for(int i=0; i<CHIP_TYPE_.size(); i++) {
            canvas->Clear();
            canvas->SetTopMargin(0.062);
            canvas->SetBottomMargin(0.14);
            canvas->SetLeftMargin(0.13);
            canvas->SetRightMargin(0.07);

            std::vector<TH1D*> hists_to_delete = {};

            TLegend* legend_main = new TLegend(0.13, 0.70, 0.43, 0.92);
            legend_main->SetFillStyle(0);
            legend_main->SetBorderSize(0);
            legend_main->SetTextSize(0.04);

            TLegend* legend_sub = new TLegend(0.33, 0.81, 0.63, 0.92);
            legend_sub->SetFillStyle(0);
            legend_sub->SetBorderSize(0);
            legend_sub->SetTextSize(0.04);

            bool is_first_plots = true;

            for(int j=0; j<CHIP_TYPE_.size(); j++) {
                if(i == j) {
                    continue;
                }
                TFile* fileSub = TFile::Open(FILENAMES_[j][3].c_str());
                TH1D* hSub = (TH1D*)fileSub->Get("h_mxAmpAC_spectra[1]");
                TH1D* hSubClone = (TH1D*)hSub->Clone(Form("hSubClone_%s_10V", CHIP_TYPE_[j].c_str()));

                for(int l=0; l<3000; l++) {
                    hSubClone->SetBinContent(l, 0);
                    hSubClone->SetBinError(l, 0);
                }

                hSubClone->Rebin(100);
                hSubClone->GetXaxis()->SetRangeUser(3000, 9000);
                if(hSubClone->GetMaximum() > 0) {
                    hSubClone->Scale(1.0 / hSubClone->GetMaximum());
                }

                // for(int l=0; l<30; l++) {
                //     hSubClone->SetBinContent(l, 0);
                //     hSubClone->SetBinError(l, 0);
                // }

                hSubClone->SetStats(0);
                hSubClone->SetTitle(";charge [ADC];normalized counts");
                hSubClone->GetYaxis()->SetRangeUser(0, 1.4);
                hSubClone->GetXaxis()->SetTitleSize(0.05); hSubClone->GetXaxis()->SetLabelSize(0.04); hSubClone->GetXaxis()->SetTitleOffset(0.8);
                hSubClone->GetYaxis()->SetTitleSize(0.05); hSubClone->GetYaxis()->SetLabelSize(0.04); hSubClone->GetYaxis()->SetTitleOffset(0.8);
                hSubClone->SetMarkerColor(kGray+2);
                hSubClone->SetLineColor(kGray+2);
                hSubClone->SetMarkerStyle(24+j);
                hSubClone->SetMarkerSize(0.8);

                hSubClone->SetLineStyle(1);
                hSubClone->SetLineWidth(2);
                if(is_first_plots) {
                    hSubClone->Draw("PE");
                    is_first_plots = false;
                } else {
                    hSubClone->Draw("SAME PE");
                }

                TF1* fSub = plot_histogram::optimise_hist_gaus(hSubClone, kGray);
                fSub->SetLineWidth(2);
                fSub->SetLineStyle(2);

                legend_sub->AddEntry(hSubClone, Form("%s, 10 V", CHIP_TYPE_[j].c_str()), "pel");
            }

            for(int k=0; k<VOLTAGE_.size(); k++) {
                TFile* fileMain = TFile::Open(FILENAMES_[i][k].c_str());
                TH1D* hMain = (TH1D*)fileMain->Get("h_mxAmpAC_spectra[1]");
                TH1D* hMainClone = (TH1D*)hMain->Clone(Form("hMainClone_%s_%sV", CHIP_TYPE_[i].c_str(), VOLTAGE_[k].c_str()));
                // threshold
                for(int l=0; l<3000; l++) {
                    hMainClone->SetBinContent(l, 0);
                    hMainClone->SetBinError(l, 0);
                }
                hMainClone->Rebin(100);
                hMainClone->GetXaxis()->SetRangeUser(3000, 9000);
                hMainClone->Scale(1.0 / hMainClone->GetMaximum());
                // if(hMainClone->GetMaximum() > 0) {
                //     hMainClone->Scale(1.0 / hMainClone->GetMaximum());
                // }
                // for(int l=0; l<30; l++) {
                //     hMainClone->SetBinContent(l, 0);
                //     hMainClone->SetBinError(l, 0);
                // }

                int color = kAzure - 3 + (2 * k);
                hMainClone->SetStats(0);
                hMainClone->SetTitle(";charge [ADC];normalized counts");
                hMainClone->GetYaxis()->SetRangeUser(0, 1.4);
                //hMainClone->GetXaxis()->SetRangeUser(3000, 9000);
                hMainClone->SetMarkerColor(color);
                hMainClone->SetLineColor(color);
                hMainClone->SetMarkerStyle(20+k);
                hMainClone->SetMarkerSize(1);
 
                hMainClone->SetLineStyle(1);
                hMainClone->SetLineWidth(2);
                if(is_first_plots) {
                    hMainClone->Draw("PE");
                    is_first_plots = false;
                } else {
                    hMainClone->Draw("SAME PE");
                }
                //hMainClone->Draw("SAME PE");

                TF1* fMain = plot_histogram::optimise_hist_gaus(hMainClone, color);
                // TF1* fMain = new TF1(Form("fMain_%s_%sV", CHIP_TYPE_[i].c_str(), VOLTAGE_[k].c_str()), "gaus", 3000, 8000);
                // hMainClone->Fit(fMain, "QNLS+","",3000, 9000);
                // fMain->SetLineStyle(1);
                // fMain->SetLineWidth(2);
                // fMain->SetRange(hMainClone->GetXaxis()->GetXmin(), hMainClone->GetXaxis()->GetXmax());
                // fMain->Draw("SAME");

                legend_main->AddEntry(hMainClone, Form("%s, %s V", CHIP_TYPE_[i].c_str(), VOLTAGE_[k].c_str()));
            }

            condition.DrawLatexNDC(0.55, 0.90, "^{55}Fe @Hiroshima Univ.(Aug. 2024)");
            condition.DrawLatexNDC(0.55, 0.86, "SeedThd 1000 ADC/cluster size = 1");
            condition.DrawLatexNDC(0.55, 0.82, Form("Plotted on %s", TIME_.c_str()));
            legend_main->Draw();
            legend_sub->Draw();
            canvas->SaveAs(Form("plot/ClusterCharge_for_labtest_%s_with_Others.pdf", CHIP_TYPE_[i].c_str()));
            //canvas->SaveAs(Form("plot/ClusterCharge_for_labtest_%s_with_Others.root", CHIP_TYPE_[i].c_str()));
        } // end of the main loop
    }

    bool charge_as_voltage = true;
    if(charge_as_voltage) {
        for(int i=0; i<VOLTAGE_.size(); i++) {
            canvas->Clear();
            canvas->SetTopMargin(0.062);
            canvas->SetBottomMargin(0.14);
            canvas->SetLeftMargin(0.13);
            canvas->SetRightMargin(0.07);

            TLegend* legend_main = new TLegend(0.13, 0.75, 0.43, 0.93);
            legend_main->SetFillStyle(0);
            legend_main->SetBorderSize(0);
            legend_main->SetTextSize(0.04);

            TLegend* legend_sub = new TLegend(0.33, 0.75, 0.63, 0.93);
            legend_sub->SetFillStyle(0);
            legend_sub->SetBorderSize(0);
            legend_sub->SetTextSize(0.04);

            bool is_first_plot = true;

            // Sub plots
            for(int j=0; j<VOLTAGE_.size(); j++) {
                if(i == j) {
                    continue;
                }
                // for(int k=0; k<CHIP_TYPE_.size(); k++) {
                //     TFile* fileSub = TFile::Open(FILENAMES_[k][j].c_str());
                //     TH1D* hSub = (TH1D*)fileSub->Get("h_mxAmpAC_spectra[1]");
                //     TH1D* hSubClone = (TH1D*)hSub->Clone(Form("hSubClone_%s_%sV", CHIP_TYPE_[k].c_str(), VOLTAGE_[j].c_str()));

                //     for(int l=0; l<3000; l++) {
                //         hSubClone->SetBinContent(l, 0);
                //         hSubClone->SetBinError(l, 0);
                //     }

                //     hSubClone->Rebin(100);
                //     hSubClone->GetXaxis()->SetRangeUser(3000, 9000);
                //     if(hSubClone->GetMaximum() > 0) {
                //         hSubClone->Scale(1.0 / hSubClone->GetMaximum());
                //     }

                //     hSubClone->SetStats(0);
                //     hSubClone->SetTitle(";charge [ADC];normalized counts");
                //     hSubClone->GetYaxis()->SetRangeUser(0, 1.4);
                //     hSubClone->GetXaxis()->SetTitleSize(0.05); hSubClone->GetXaxis()->SetLabelSize(0.04); hSubClone->GetXaxis()->SetTitleOffset(0.8);
                //     hSubClone->GetYaxis()->SetTitleSize(0.05); hSubClone->GetYaxis()->SetLabelSize(0.04); hSubClone->GetYaxis()->SetTitleOffset(0.8);            
                //     hSubClone->SetMarkerColor(kGray+2);
                //     hSubClone->SetLineColor(kGray+2);
                //     hSubClone->SetMarkerStyle(24+j);
                //     hSubClone->SetMarkerSize(0.8);
                //     hSubClone->SetLineStyle(1);
                //     hSubClone->SetLineWidth(2);

                //     if(is_first_plot) {
                //         hSubClone->Draw("PE");
                //         is_first_plot = false;
                //     } else {
                //         hSubClone->Draw("SAME PE");
                //     }

                //     TF1* fSub = plot_histogram::optimise_hist_gaus(hSubClone, kGray);
                //     fSub->SetLineWidth(2);
                //     fSub->SetLineStyle(2);

                //     legend_sub->AddEntry(hSubClone, Form("%s, %s V", CHIP_TYPE_[k].c_str(), VOLTAGE_[j].c_str()), "pel");
                // }

                TFile* fileSub = TFile::Open(FILENAMES_[1][j].c_str());
                TH1D* hSub = (TH1D*)fileSub->Get("h_mxAmpAC_spectra[1]");
                TH1D* hSubClone = (TH1D*)hSub->Clone(Form("hSubClone_std_%sV", VOLTAGE_[j].c_str()));

                for(int l=0; l<3000; l++) {
                    hSubClone->SetBinContent(l, 0);
                    hSubClone->SetBinError(l, 0);
                }

                hSubClone->Rebin(100);
                hSubClone->GetXaxis()->SetRangeUser(3000, 9000);
                if(hSubClone->GetMaximum() > 0) {
                    hSubClone->Scale(1.0 / hSubClone->GetMaximum());
                }

                hSubClone->SetStats(0);
                hSubClone->SetTitle(";charge [ADC];normalized counts");
                hSubClone->GetYaxis()->SetRangeUser(0, 1.4);
                hSubClone->GetXaxis()->SetTitleSize(0.05); hSubClone->GetXaxis()->SetLabelSize(0.04); hSubClone->GetXaxis()->SetTitleOffset(0.8);
                hSubClone->GetYaxis()->SetTitleSize(0.05); hSubClone->GetYaxis()->SetLabelSize(0.04); hSubClone->GetYaxis()->SetTitleOffset(0.8);            
                hSubClone->SetMarkerColor(kGray+2);
                hSubClone->SetLineColor(kGray+2);
                hSubClone->SetMarkerStyle(24+j);
                hSubClone->SetMarkerSize(0.8);
                hSubClone->SetLineStyle(1);
                hSubClone->SetLineWidth(2);

                if(is_first_plot) {
                    hSubClone->Draw("PE");
                    is_first_plot = false;
                } else {
                    hSubClone->Draw("SAME PE");
                }

                TF1* fSub = plot_histogram::optimise_hist_gaus(hSubClone, kGray);
                fSub->SetLineWidth(2);
                fSub->SetLineStyle(2);
                legend_sub->AddEntry(hSubClone, Form("std, %s V", VOLTAGE_[j].c_str()), "pel");
            } // VOLATGE_

            // Main plots
            for(int k=0; k<CHIP_TYPE_.size(); k++) {
                TFile* fileMain = TFile::Open(FILENAMES_[k][i].c_str());
                TH1D* hMain = (TH1D*)fileMain->Get("h_mxAmpAC_spectra[1]");
                TH1D* hMainClone = (TH1D*)hMain->Clone(Form("hMainClone_%s_%sV", CHIP_TYPE_[k].c_str(), VOLTAGE_[i].c_str()));

                for(int l=0; l<3000; l++) {
                    hMainClone->SetBinContent(l, 0);
                    hMainClone->SetBinError(l, 0);
                }

                hMainClone->Rebin(100);
                hMainClone->GetXaxis()->SetRangeUser(3000, 9000);
                if(hMainClone->GetMaximum() > 0) {
                    hMainClone->Scale(1.0 / hMainClone->GetMaximum());
                }

                int color = kAzure - 3 + (2 * k);
                hMainClone->SetStats(0);
                hMainClone->SetTitle(";charge [ADC];normalized counts");
                hMainClone->GetYaxis()->SetRangeUser(0, 1.4);
                hMainClone->GetXaxis()->SetTitleSize(0.05); hMainClone->GetXaxis()->SetLabelSize(0.04); hMainClone->GetXaxis()->SetTitleOffset(0.8);
                hMainClone->GetYaxis()->SetTitleSize(0.05); hMainClone->GetYaxis()->SetLabelSize(0.04); hMainClone->GetYaxis()->SetTitleOffset(0.8);            
                hMainClone->SetMarkerColor(color);
                hMainClone->SetLineColor(color);
                hMainClone->SetMarkerStyle(24+k);
                hMainClone->SetMarkerSize(1);

                hMainClone->SetLineStyle(1);
                hMainClone->SetLineWidth(2);
                if(is_first_plot) {
                    hMainClone->Draw("PE");
                } else {
                    hMainClone->Draw("SAME PE");
                }

                TF1* fMain = plot_histogram::optimise_hist_gaus(hMainClone, color);

                legend_main->AddEntry(hMainClone, Form("%s, %s V", CHIP_TYPE_[k].c_str(), VOLTAGE_[i].c_str()), "pel");
            } // CHIP_TYPE_

            condition.DrawLatexNDC(0.55, 0.90, "^{55}Fe @Hiroshima Univ.(Aug. 2024)");
            condition.DrawLatexNDC(0.55, 0.86, "SeedThd 1000 ADC/cluster size = 1");
            condition.DrawLatexNDC(0.55, 0.82, Form("Plotted on %s", TIME_.c_str()));
            legend_main->Draw();
            legend_sub->Draw();
            canvas->SaveAs(Form("plot/ClusterCharge_for_labtest_%sV_with_Others.pdf", VOLTAGE_[i].c_str()));
            //canvas->SaveAs(Form("plot/ClusterCharge_for_labtest_%sV_with_Others.root", VOLTAGE_[i].c_str()));
        } // end of the main loop
    }

    bool plot_voltage_dependence = true;
    if(plot_voltage_dependence) {
        std::vector<std::vector<double>> vSigma;
        std::vector<std::vector<double>> vMean;
        std::vector<std::vector<double>> vSigmaError;
        std::vector<std::vector<double>> vMeanError;
        
        vSigma.resize(CHIP_TYPE_.size());
        vMean.resize(CHIP_TYPE_.size());
        vSigmaError.resize(CHIP_TYPE_.size());
        vMeanError.resize(CHIP_TYPE_.size());
        for(size_t int i=0; i<VOLTAGE_.size(); i++) {
            vSigma[i].resize(VOLTAGE_.size());
            vMean[i].resize(VOLTAGE_.size());
            vSigmaError[i].resize(VOLTAGE_.size());
            vMeanError[i].resize(VOLTAGE_.size());
        }

        // Sigma and Mean
        for(int i=0; i<CHIP_TYPE_.size(); i++) {
            for(int j=0; j<VOLTAGE_.size(); j++) {
                TFile* fileMain = TFile::Open(FILENAMES_[i][j].c_str());
                TH1D* hMain = (TH1D*)fileMain->Get("h_mxAmpAC_spectra[1]");
                TH1D* hMainClone = (TH1D*)hMain->Clone(Form("hMainClone_%s_%sV", CHIP_TYPE_[i].c_str(). VOLTAGE_[k].c_str()));

                for(int l=0; l<3000; l++) {
                    hMainClone->SetBinContent(l, 0);
                    hMainClone->SetBinError(l, 0);
                }

                TF1* fMain = plot_histogram::optimise_hist_gaus(hMainClone, 1);
                vSigma[i][j].push_back(fMain->GetParameter(2));
                vMean[i][j].push_back(fMain->GetParameter(1));
                vSigmaError[i][j].push_back(fMain->GetParError(2));
                vMeanError[i][j].push_back(fMain->GetParError(1));
            } // VOLTAGE_
        } // CHIP_TYPE_

        // Sigma
        canvas->Clear();
        canvas->SetTopMargin(0.062);
        canvas->SetBottomMargin(0.14);
        canvas->SetLeftMargin(0.13);
        canvas->SetRightMargin(0.07);

        TLegend* legend_main = new TLegend(0.13, 0.70, 0.43, 0.92);
        legend_main->SetFillStyle(0);
        legend_main->SetBorderSize(0);
        legend_main->SetTextSize(0.04);

        TLegend* legend_sub = new TLegend(0.33, 0.81, 0.63, 0.92);
        legend_sub->SetFillStyle(0);
        legend_sub->SetBorderSize(0);
        legend_sub->SetTextSize(0.04);

        

    }
}