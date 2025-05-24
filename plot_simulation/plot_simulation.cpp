#ifndef PLOT_SIMULATION_H
#define PLOT_SIMULATION_H

#include "plot_simulation.h"

// set histogram style for ROOT
void plot_simulation::set_rootStyle() {
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "---------------------------- Set ROOT style -----------------------------" << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTopMargin(0.12); // changed from 0.02
    gStyle->SetPadBottomMargin(0.14); // changed from 0.12
    gStyle->SetPadLeftMargin(0.12); // changed from 0.14
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetMarkerColor(kRed);
    gStyle->SetMarkerSize(1.5);
    gStyle->SetMarkerStyle(20);
    gStyle->SetLineWidth(2);
    gStyle->SetFuncColor(kBlack);
    gStyle->SetFuncWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.045,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.06,"xyz");
    gStyle->SetTitleOffset(0.00,"y");
    gStyle->SetTitleOffset(0.8,"x");
    gStyle->SetTitleOffset(0.00,"z");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    gStyle->SetOptStat(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.04);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptFit(0);
}

// set cluster size distribution style
void plot_simulation::set_clSizeStyle(TH1D* hist, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle("cluster size");
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(0, range);
}

// set correlation distribution style: ref vs hit
void plot_simulation::set_correlationStyle(TH1D* hist, const std::string& axis, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle(Form("%s_{ref}-%s[mm]", axis.c_str(), axis.c_str()));
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(-range, range);
}

// set charge distribution style
void plot_simulation::set_chargeStyle(TH1D* hist, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle("charge [ADC]");
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(0, range);
}

// set residual distribution style
void plot_simulation::set_residualStyle(TH1D* hist, const std::string& axis, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle(Form("%s_{ref}-%s[mm]", axis.c_str(), axis.c_str()));
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(-range, range);
}

// set track chi2 distribution style
void plot_simulation::set_trackChi2Style(TH1D* hist, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle("#chi^{2}");
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(0, range);
}

// set cluster multiplicity per event distribution style
void plot_simulation::set_multiplicityStyle(TH1D* hist, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle("clusters");
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(0, range);
}

// set track chi2/ndof distribution style
void plot_simulation::set_trackChi2ndofStyle(TH1D* hist, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle("#chi^{2}/NDF");
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(0, range);
}

// set track angle distribution style
void plot_simulation::set_trackAngleStyle(TH1D* hist, const std::string& axis, double range, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle(Form("angle_%s [rad]", axis.c_str()));
    hist->GetYaxis()->SetTitle("events");
    hist->GetXaxis()->SetRangeUser(-range, range);
}

// set hitmap projection distribution style
void plot_simulation::set_projectionStyle(TH1D* hist, const std::string& axis, int rebin) {
    hist->Rebin(rebin);
    hist->GetXaxis()->SetTitle(Form("%s[px]", axis.c_str()));
    hist->GetYaxis()->SetTitle("events");
}

// set TH1D histogram style
// change histogram style based on the histogram title
void plot_simulation::set_th1dStyle(TH1D* hist, int color) {
    std::string hist_title = hist->GetTitle();
    std::cout << "Set TH1D style(plot_simulation) -" << hist->GetName() << std::endl;
    std::cout << "Hist title: " << hist_title << std::endl;

    // set histogram global style
    hist->Sumw2();
    hist->SetStats(0);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.5);

    // set title and label seze and position
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleOffset(0.8);

    // set title and range based on the histogram title
    //TODO: ALPIDE residuals units should be converted to um
    //TODO: Add scaled histogram GetMaximum() or Integral()
    // ALPIDE histogram style
    if(search_word(hist_title, "ALPIDE")) {
        // std::cout << "ALPIDE histogram is changed" << std::endl;
        // ALPIDE: cluster size
        if(search_word(hist_title, "cluster size")) {
            set_clSizeStyle(hist, 10, 1);
        }
        //ALPIDE: correlation
        if(search_word(hist_title, "correlation")) {
            if(search_word(hist_title, "X")) {
                set_correlationStyle(hist, "x", 6, 10);
            }
            else if(search_word(hist_title, "Y")) {
                set_correlationStyle(hist, "y", 6, 10);
            }
        }
        // ALPIDE: charge distribution ex. seed charge, cluster charge, neighbor charge
        if(search_word(hist_title, "charge")) {
            set_chargeStyle(hist, 5000, 1);
        }
        // ALPIDE: residual distribution
        if(search_word(hist_title, "residual")) {
            if(search_word(hist_title, "X")) {
                set_residualStyle(hist, "x", 60, 10);
            }
            else if(search_word(hist_title, "Y")) {
                set_residualStyle(hist, "y", 60, 10);
            }
        }
        // ALPIDE: cluster multiplicity per event
        if(search_word(hist_title, "multiplicity")) {
            set_multiplicityStyle(hist, 10, 1);
        }
        // ALPIDE: hitmap projection
        if(search_word(hist_title, "position")) {
            if(search_word(hist->GetXaxis()->GetTitle(), "x")) {
                set_projectionStyle(hist, "x", 10);
            }
            else if(search_word(hist->GetXaxis()->GetTitle(), "y")) {
                set_projectionStyle(hist, "y", 10);
            }
        }
    }
    // CE65 histogram style
    if(search_word(hist_title, "CE65")) {
        // std::cout << "CE65 histogram style is changed" << std::endl;
        // CE65: cluster size
        if(search_word(hist_title, "cluster size")) {
            set_clSizeStyle(hist, 100, 1);
        }
        // CE65: correlation
        if(search_word(hist_title, "correlation")) {
            if(search_word(hist_title, "X")) {
                set_correlationStyle(hist, "x", 6, 10);
            }
            else if(search_word(hist_title, "Y")) {
                set_correlationStyle(hist, "y", 6, 10);
            }
        }
        // CE65: charge distribution ex. seed charge, cluster charge, neighbor charge
        if(search_word(hist_title, "charge")) {
            set_chargeStyle(hist, 5000, 10);
        }
        // CE65: hitmap projection
        if(search_word(hist_title, "position")) {
            if(search_word(hist->GetXaxis()->GetTitle(), "x")) {
                set_projectionStyle(hist, "x", 10);
            }
            else if(search_word(hist->GetXaxis()->GetTitle(), "y")) {
                set_projectionStyle(hist, "y", 10);
            }
        }
    }

    // CE65: residual distribution
    if(search_word(hist_title, "Residual in local")) {
            if(search_word(hist_title, "X")) {
                set_residualStyle(hist, "x", 60, 10);
            }
            else if(search_word(hist_title, "Y")) {
                set_residualStyle(hist, "y", 60, 10);
            }
        }

    // Track histogram style
    if(search_word(hist_title, "Track")) {
        // Track: chi2/ndof
        if(search_word(hist_title, "chi2ndof")) {
            set_trackChi2ndofStyle(hist, 10, 1);
        }
        // Track: chi2
        if(search_word(hist_title, "chi2")) {
            set_trackChi2Style(hist, 100, 1);
        }
        // Track: track angle
        if(search_word(hist_title, "angle")) {
            if(search_word(hist_title, "X")) {
                set_trackAngleStyle(hist, "x", 0.1, 20);
            }
            else if(search_word(hist_title, "Y")) {
                set_trackAngleStyle(hist, "y", 0.1, 20);
            }
        }
    }
}

// TH2D histogram style: hitmap
void plot_simulation::set_th2dStyle(TH2D* hist) {
    std::cout << "Set TH2D style -" << hist->GetName() << std::endl;
    hist->SetStats(0);
}

TH1D* plot_simulation::search_TH1D_inCanvas(TCanvas* canvas) {
    // get list of primitives in canvas
    if(canvas == nullptr) {
        std::cout << "Canvas is not found" << std::endl;
        return nullptr;
    }
    TIter next(canvas->GetListOfPrimitives());
    TObject* obj = nullptr;
    while ((obj = next())) {
        // search TH1 object in canvas
        if (obj->InheritsFrom(TH1D::Class())) {
            TH1D* hist = dynamic_cast<TH1D*>(obj);
            std::cout << "TH1D is found -" << obj->GetName() << std::endl;
            return hist;
        }
        // // else if(obj->InheritsFrom(TPad::Class())) {
        // //     TPad* pad = dynamic_cast<TPad*>(obj);
        // //     TH1D* hist = search_TH1D_inCanvas(pad);
        // //     if(hist != nullptr) {
        // //         return hist;
        // //     }
        // // }
    }
    std::cout << "TH1D is not found" << std::endl;
    return nullptr;
}

// search TF1 object in canvas or TH1 object
TF1* plot_simulation::search_TF1_inCanvas(TCanvas* canvas) {
    // get list of primitives in canvas
    TIter next(canvas->GetListOfPrimitives());
    TObject* obj = nullptr;
    while ((obj = next())) {
        // search TF1 object in canvas
        if (obj->InheritsFrom(TF1::Class())) {
            TF1* fit = dynamic_cast<TF1*>(obj);
            std::cout << "TF1 is found -" << obj->GetName() << std::endl;
            return fit;
        }
        // search TH1 object in canvas
        else
        if (obj->InheritsFrom(TH1::Class())) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            std::cout << "TH1 is found -" << obj->GetName() << std::endl;
            TList* obj_list = hist->GetListOfFunctions();
            // search TF1 object in TH1 object
            if(obj_list) {
                TIter funcIter(obj_list);
                TObject* funcObj = nullptr;
                while ((funcObj = funcIter())) {
                    if (funcObj->InheritsFrom(TF1::Class())) {
                        TF1* fit = dynamic_cast<TF1*>(funcObj);
                        std::cout << "TF1 is found -" << funcObj->GetName() << std::endl;
                        return fit;
                    }
                }
            }
        }
    }
    std::cout << "TF1 is not found" << std::endl;
    return nullptr;
}

void plot_simulation::write_plot(std::variant<TH1D*, TH2D*> hist, const std::string& hist_name) {
    std::cout << "Write plot to canvas : " << hist_name << std::endl;
    
    std::visit([&hist_name](auto* hist) {
        TCanvas *c = new TCanvas(hist->GetName(), hist->GetName(), 1920, 1080);
        if (!hist) {
            print_message("Histogram is not found", RED);
            return;
        }
        if constexpr (std::is_same_v<std::decay_t<decltype(hist)>, TH1D*>) {
            std::cout << "Histogram type: TH1D" << std::endl;
            hist->Draw("PE");
            if(search_word(hist->GetTitle(), "clustercharge")) {
                if(search_TF1_inCanvas(c)!=nullptr) {
                    draw_langausFit_result(hist, search_TF1_inCanvas(c));
                }
            }
            else {
                if(search_TF1_inCanvas(c)!=nullptr) {
                    draw_gausFit_result(hist, search_TF1_inCanvas(c));
                }
                else {
                    draw_th1d_result(hist);
                }
            }
        } else if constexpr (std::is_same_v<std::decay_t<decltype(hist)>, TH2D*>) {
            std::cout << "Histogram type: TH2D" << std::endl;
            hist->Draw("colz");
        }
        c->Write(hist_name.c_str());
    }, hist);
}

// for 6ALPIDE
// TODO: need to be updated for refname and dutname
// TODO: to touch TCanvas and search TH1D object
// This code cannot find some item
void plot_simulation::write_allRefProjection(TFile* input_file, TFile* output_file, std::vector<std::string> refnames, std::string axis) {
    std::cout << "Write all reference projection" << std::endl;

    std::string input_dirName, output_dirName;
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1920, 1080);

    input_file->ls();
    TIter next(input_file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        std::cout << "Key name: " << key->GetName() << std::endl;
    }

    for(int i=0; i<refnames.size(); i++) {
        input_dirName = "Hitmap/" + refnames[i] + "/";
        output_dirName = "Hitmap/";

        // input_file->cd(input_dirName.c_str());
        // input_file->ls();

        // TH1D* hProjection = (TH1D*)input_file->Get((input_dirName + "hProjection" + axis).c_str());
        std::cout << "Get histogram from " << input_dirName + "hProjection" + axis << std::endl;
        TCanvas *canvas_histogram = (TCanvas*)input_file->Get((input_dirName + "hProjection" + axis).c_str());
        TH1D* hProjection = search_TH1D_inCanvas(canvas_histogram);

        if(hProjection == nullptr) {
            print_message(Form("Histogram is not found at %s", input_dirName.c_str()), RED);
            return;
        }
        set_projectionStyle(hProjection, axis, 10);

        if(i == 0) {
            std::cout << "Draw first histogram :" << hProjection->GetName() << std::endl;
            hProjection->Draw("PE");
        }
        else {
            std::cout << "Draw next histogram :" << hProjection->GetName() << std::endl;
            hProjection->Draw("PEsame");
        }
    }
    // canvas->Write(("hProjection" + axis).c_str());
    output_file->cd(output_dirName.c_str());
    canvas->Write(("hProjection" + axis).c_str());
    
}

void plot_simulation::run() {
    set_rootStyle();

    refnames = {"ALPIDE_0", "ALPIDE_1", "ALPIDE_2", "ALPIDE_4", "ALPIDE_5", "ALPIDE_6"};
    dutnames = {"CE65_3"};
    input_fileDir = "/home/towa/alice3/hist/kek_sim/";
    // base_filename = "05_analysis_kek_sim_CE65_p15_gap_6ALPIDE_202502131501.root";
    // base_filename = "05_analysis_kek_sim_CE65_p15_std_6ALPIDE_202502131448.root";
    // base_filename = "05_analysis_kek_sim_CE65_p22p5_gap_6ALPIDE_202502131529.root";
    base_filename = "05_analysis_kek_sim_CE65_p22p5_std_6ALPIDE_202502131515.root";
    input_filename = input_fileDir + base_filename;
    output_filename = "plot_" + base_filename;
    rootFile_dir_names = {"Hitmap", "Correlation", "AnalysisRef", "Track", "AnalysisCE65"};

    input_file = new TFile(input_filename.c_str());
    output_file = new TFile(output_filename.c_str(), "RECREATE");
    TFile *print_file = new TFile(("print_" + output_filename).c_str(), "RECREATE");

    hitmap = output_file->mkdir("Hitmap");
    correlation = output_file->mkdir("Correlation");
    analysisRef = output_file->mkdir("AnalysisRef");
    track = output_file->mkdir("Track");
    analysisCE65 = output_file->mkdir("AnalysisCE65");

    for (const auto refname : refnames) {
        TDirectory *ref_dir = hitmap->mkdir(refname.c_str());
        TDirectory *ref_dir_correlation = correlation->mkdir(refname.c_str());
        TDirectory *ref_dir_analysis = analysisRef->mkdir(refname.c_str());
        plot_histogram::draw_ref_clusteringSpatial(refname, input_file, output_file);
        plot_histogram::draw_correlation(refname, input_file, output_file);
        plot_histogram::draw_ref_analysis(refname, input_file, output_file);
    }

    for(const auto dutname : dutnames) {
        TDirectory *dut_dir = hitmap->mkdir(dutname.c_str());
        TDirectory *dut_dir_correlation = correlation->mkdir(dutname.c_str());
        TDirectory *dut_dir_analysis = analysisCE65->mkdir(dutname.c_str());
        plot_histogram::draw_dut_clusteringAnalog(dutname, input_file, output_file);
        plot_histogram::draw_correlation(dutname, input_file, output_file);
        plot_histogram::draw_AnalysisCE65(dutname, input_file, output_file);
    }

    plot_histogram::draw_Tracking4D(input_file, output_file);
    output_file->Close();

    // std::cout << "File name: " << output_file->GetName() << std::endl;
    // plot_simulation::write_allRefProjection(output_file, print_file, refnames, "X");
    // std::cout << "Finished plot_simulation::run()" << std::endl;
}

plot_simulation::plot_simulation() {
    LOG_DEBUG.source("plot_simulation::plot_simulation()") << "Initialize plot_simulation.";
}

#endif // PLOT_SIMULATION_H