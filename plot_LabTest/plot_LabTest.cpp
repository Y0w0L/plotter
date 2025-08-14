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
                hMainClone->SetMarkerStyle(20+k);
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
        for(int i=0; i<CHIP_TYPE_.size(); i++) {
            vSigma[i].resize(VOLTAGE_.size());
            vMean[i].resize(VOLTAGE_.size());
            vSigmaError[i].resize(VOLTAGE_.size());
            vMeanError[i].resize(VOLTAGE_.size());
        }

        std::vector<double> vVoltage;
        std::vector<double> vVoltageError;
        for(size_t i=0; i<VOLTAGE_.size(); i++) {
            vVoltage.push_back(std::stod(VOLTAGE_[i]));
            vVoltageError.push_back(0.0);
        }

        // Sigma and Mean
        for(int i=0; i<CHIP_TYPE_.size(); i++) {
            for(int j=0; j<VOLTAGE_.size(); j++) {
                TFile* fileMain = TFile::Open(FILENAMES_[i][j].c_str());
                TH1D* hMain = (TH1D*)fileMain->Get("h_mxAmpAC_spectra[1]");
                TH1D* hMainClone = (TH1D*)hMain->Clone(Form("hMainClone_%s_%sV", CHIP_TYPE_[i].c_str(), VOLTAGE_[j].c_str()));

                for(int l=0; l<3000; l++) {
                    hMainClone->SetBinContent(l, 0);
                    hMainClone->SetBinError(l, 0);
                }

                TF1* fMain = plot_histogram::optimise_hist_gaus(hMainClone, 1);
                // vSigma[i].push_back(fMain->GetParameter(2));
                // vMean[i].push_back(fMain->GetParameter(1));
                // vSigmaError[i].push_back(fMain->GetParError(2));
                // vMeanError[i].push_back(fMain->GetParError(1));

                vSigma[i][j] = fMain->GetParameter(2);
                vSigmaError[i][j] = fMain->GetParError(2);
                vMean[i][j] = fMain->GetParameter(1);
                vMeanError[i][j] = fMain->GetParError(1);

            } // VOLTAGE_
        } // CHIP_TYPE_

        // Sigma
        canvas->Clear();
        canvas->SetTopMargin(0.062);
        canvas->SetBottomMargin(0.14);
        canvas->SetLeftMargin(0.13);
        canvas->SetRightMargin(0.07);

        TLegend* legend_sigma = new TLegend(0.8, 0.70, 0.92, 0.92);
        legend_sigma->SetFillStyle(0);
        legend_sigma->SetBorderSize(0);
        legend_sigma->SetTextSize(0.04);

        TLegend* legend_mean = new TLegend(0.8, 0.50, 0.92, 0.82);
        legend_mean->SetFillStyle(0);
        legend_mean->SetBorderSize(0);
        legend_mean->SetTextSize(0.04);

        // plot for sigma vs voltage
        TMultiGraph* mg_sigma = new TMultiGraph();
        mg_sigma->SetTitle(";voltage [V];peak width");

        int colors[] = {kRed+2, kAzure+2, kGreen+2};
        int line_colors[] = {kRed-9, kAzure-9, kGreen-9};
        int markers[] = {20,21,22};

        for(int i=0; i<CHIP_TYPE_.size(); i++) {
            TGraphErrors* gr_sigma = new TGraphErrors(vVoltage.size(), &vVoltage[0], &vSigma[i][0], &vVoltageError[0], &vSigmaError[i][0]);
            gr_sigma->SetMarkerColor(colors[i]);
            gr_sigma->SetLineColor(colors[i]);
            gr_sigma->SetMarkerStyle(markers[i]);
            //gr_sigma->SetLineStyle(i+1);
            //gr_sigma->SetLineWidth(0);
            gr_sigma->SetMarkerSize(1.5);
            //mg_sigma->Add(gr_sigma, "EP");

            TGraphErrors* gr_sigma_line = new TGraphErrors(vVoltage.size(), &vVoltage[0], &vSigma[i][0], &vVoltageError[0], &vSigmaError[i][0]);
            gr_sigma_line->SetLineColor(line_colors[i]);
            //gr_sigma_line->SetLineColor(kGray);
            gr_sigma_line->SetLineStyle(2);
            gr_sigma_line->SetLineWidth(2);

            mg_sigma->Add(gr_sigma_line, "L");
            mg_sigma->Add(gr_sigma, "EP");

            legend_sigma->AddEntry(gr_sigma, Form("%s", CHIP_TYPE_[i].c_str()), "pe");
        }

        mg_sigma->Draw("A");
        mg_sigma->GetXaxis()->SetTitleSize(0.05); mg_sigma->GetXaxis()->SetLabelSize(0.04);// mg_sigma->SetTitleOffset(0.8);
        mg_sigma->GetYaxis()->SetTitleSize(0.05); mg_sigma->GetYaxis()->SetLabelSize(0.04);// mg_sigma->SetTitleOffset(0.8);
        //mg_sigma->Draw();

        legend_sigma->Draw();
        condition.DrawLatexNDC(0.58, 0.25, "^{55}Fe @Hiroshima Univ.(Aug. 2024)");
        condition.DrawLatexNDC(0.58, 0.21, "SeedThd 1000 ADC/cluster size = 1");
        condition.DrawLatexNDC(0.58, 0.17, Form("Plotted on %s", TIME_.c_str()));

        canvas->SaveAs("plot/Sigma_vs_Voltage_as_labtest.pdf");

        canvas->Clear();
        // plot for mean vs voltage
        TMultiGraph* mg_mean = new TMultiGraph();
        mg_mean->SetTitle(";voltage [V];peak position");
        
        for(int i=0; i<CHIP_TYPE_.size(); i++) {
            TGraphErrors* gr_mean = new TGraphErrors(vVoltage.size(), &vVoltage[0], &vMean[i][0], &vVoltageError[0], &vMeanError[i][0]);
            TGraphErrors* gr_mean_line = new TGraphErrors(vVoltage.size(), &vVoltage[0], &vMean[i][0], &vVoltageError[0], &vMeanError[i][0]);
            
            gr_mean->SetMarkerColor(colors[i]);
            gr_mean->SetLineColor(colors[i]);
            gr_mean->SetMarkerStyle(markers[i]);
            //gr_sigma->SetLineStyle(i+1);
            //gr_sigma->SetLineWidth(0);
            gr_mean->SetMarkerSize(1.5);
            
            gr_mean_line->SetLineColor(line_colors[i]);
            gr_mean_line->SetLineStyle(2);
            gr_mean_line->SetLineWidth(2);

            mg_mean->Add(gr_mean_line, "L");
            mg_mean->Add(gr_mean, "EP");

            legend_mean->AddEntry(gr_mean, Form("%s", CHIP_TYPE_[i].c_str()), "pe");
        }

        mg_mean->Draw("A");
        mg_mean->GetXaxis()->SetTitleSize(0.05); mg_mean->GetXaxis()->SetLabelSize(0.04);// mg_sigma->SetTitleOffset(0.8);
        mg_mean->GetYaxis()->SetTitleSize(0.05); mg_mean->GetYaxis()->SetLabelSize(0.04);// mg_sigma->SetTitleOffset(0.8);

        legend_mean->Draw();
        condition.DrawLatexNDC(0.58, 0.25, "^{55}Fe @Hiroshima Univ.(Aug. 2024)");
        condition.DrawLatexNDC(0.58, 0.21, "SeedThd 1000 ADC/cluster size = 1");
        condition.DrawLatexNDC(0.58, 0.17, Form("Plotted on %s", TIME_.c_str()));

        canvas->SaveAs("plot/Mean_vs_Voltage_as_labtest.pdf");
    }
}