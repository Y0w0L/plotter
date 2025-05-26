#include "plot_histogram.h"
#include "tools/langaus.C"
#include "tools/twogaus.C"

// // constructor
// plot_histogram::plot_histogram() {
//     std::cout << "plot_histogram object is created" << std::endl;
// }

plot_histogram::plot_histogram()
    : hitmap(nullptr), correlation(nullptr), analysisRef(nullptr), track(nullptr), analysisCE65(nullptr),
      input_file(nullptr), output_file(nullptr),
      hHitmapLocal(nullptr), hHitmapGlobal(nullptr), hClusterPositionLocal(nullptr),
      hProjectionX(nullptr), hProjectionY(nullptr),
      hCorrelationX(nullptr), hCorrelationY(nullptr), hCorrelationX_2D(nullptr), hCorrelationY_2D(nullptr),
      hTrackChi2(nullptr), hTrackChi2ndof(nullptr), hTrackAngleX(nullptr), hTrackAngleY(nullptr),
      hClusterSize(nullptr), hClusterSize_SeedCut(nullptr), hClusterSeedCharge(nullptr), hClusterCharge(nullptr), hClusterNeighborCharge(nullptr), hClusterNeighborChargeSum(nullptr),
      hCutHisto(nullptr), hClusterMultiplicity(nullptr),
      hResidualsX(nullptr), hResidualsY(nullptr),
      fProjectionX(nullptr), fProjectionY(nullptr),
      fCorrelationX(nullptr), fCorrelationY(nullptr),
      fTrackAngleX(nullptr), fTrackAngleY(nullptr),
      fResidualsX(nullptr), fResidualsY(nullptr), fClusterCharge(nullptr)
{
    // constructor
}
plot_histogram::~plot_histogram() {
    // destructor
    LOG_INFO << "plot_histogram object destroyed.";
}

// set histogram style for ROOT
void plot_histogram::set_rootStyle() {
    LOG_INFO.source("plot_histogram::set_rootStyle") << "-------------------------------------------------------------------------";
    LOG_INFO.source("plot_histogram::set_rootStyle") << "---------------------------- Set ROOT style -----------------------------";
    LOG_INFO.source("plot_histogram::set_rootStyle") << "-------------------------------------------------------------------------";
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
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
    gStyle->SetTitleOffset(0.6,"y");
    gStyle->SetTitleOffset(1.0,"x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    gStyle->SetOptStat(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
}

// set TH1D histogram style
void plot_histogram::set_th1dStyle(TH1D* hist, int color) {
    LOG_DEBUG.source("plot_histogram::set_th1dStyle") << "Set TH1D style -" << hist->GetName();
    hist->SetStats(0);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.5);
}

// set TH2D histogram style
void plot_histogram::set_th2dStyle(TH2D* hist) {
    LOG_DEBUG.source("plot_histogram::set_th2dStyle") << "Set TH2D style -" << hist->GetName();
    hist->SetStats(0);
}

// set TF1 fit function style
void plot_histogram::set_tf1Style(TF1* fit, int color) {
    LOG_DEBUG.source("plot_histogram::set_tf1Style") << "Set TF1 style -" << fit->GetName();
    fit->SetLineColor(color);
}

void plot_histogram::set_tlatexStyle(TLatex* latex) {
    latex->SetNDC(true);
    latex->SetTextSize(0.04);
    latex->SetTextFont(42);
    latex->SetTextColor(kBlack);
    latex->SetTextAlign(13);
}

void plot_histogram::remove_whiteSpace(std::string& word) {
    word.erase(std::remove_if(word.begin(), word.end(), ::isspace), word.end());
}

void plot_histogram::replace_string(std::string& word, const std::string& target, const std::string& replacement) {
    size_t pos = 0;
    while ((pos = word.find(target, pos)) != std::string::npos) {
        word.replace(pos, target.length(), replacement);
        pos += target.length();
    }
}

// search word in string word is search word, target is search target
bool plot_histogram::search_word(const std::string& word, const std::string& target) {
    std::string word_lower = word;
    std::string target_lower = target;
    
    remove_whiteSpace(word_lower);
    remove_whiteSpace(target_lower);

    transform(word_lower.begin(), word_lower.end(), word_lower.begin(), ::tolower);
    transform(target_lower.begin(), target_lower.end(), target_lower.begin(), ::tolower);
    
    return word_lower.find(target_lower) != std::string::npos;
}

// print message with color
void plot_histogram::print_message(std::string message, const char* color) {
    std::cout << color << message << RESET << std::endl;
}

std::string plot_histogram::format_uncertainty(double value, double error) {
    if (error == 0) {
        return std::to_string(value) + " #pm 0";
    }

    int error_exponent = static_cast<int>(std::floor(std::log10(std::abs(error))));
    double first_digit = error / std::pow(10, error_exponent);
    double rounded_error = std::round(first_digit) * std::pow(10, error_exponent);

    double rounded_value = std::round(value / std::pow(10, error_exponent)) * std::pow(10, error_exponent);

    int common_exponent = std::min(error_exponent, static_cast<int>(std::floor(std::log10(std::abs(value)))));
    
    if (common_exponent < -2 || common_exponent > 3) {
        double scaled_value = rounded_value / std::pow(10, common_exponent);
        double scaled_error = rounded_error / std::pow(10, common_exponent);

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << scaled_value
            << " #pm " << std::fixed << std::setprecision(3) << scaled_error;
        return oss.str();
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(0) << rounded_value
        << " #pm " << rounded_error;
    return oss.str();
}

// write histogram to ROOT file
void plot_histogram::write_plot(std::variant<TH1D*, TH2D*> hist, const std::string& hist_name) {
    LOG_INFO.source("plot_histogram::write_plot") << "Write plot : " << hist_name;
    std::visit([&hist_name](auto* hist) {
        if(!hist) {
            LOG_ERROR.source("plot_histogram::write_plot") << "Histogram is not found";
            throw std::invalid_argument("Input histogram class is empty.");
        }
        hist->Write(hist_name.c_str());
    }, hist);
}

std::string plot_histogram::get_tf1Parameter(TF1* fit, int par) {
    std::string par_name;
    std::string par_value;
    double value, value_error;

    par_name = fit->GetParName(par);
    value = fit->GetParameter(par);
    value_error = fit->GetParError(par);

    par_value = format_uncertainty(value, value_error);
    //return par_name + " = " + std::to_string(value) + " #pm " + std::to_string(value_error);
    return Form("%s = %s", par_name.c_str(), par_value.c_str());
}

std::string plot_histogram::get_chi2ndfParameter(TF1* fit) {
    double chi2, ndf, chi2_ndf;

    chi2 = fit->GetChisquare();
    ndf = fit->GetNDF();
    chi2_ndf = chi2 / ndf;

    //return "#chi^{2}/NDF = " + std::to_string(chi2) + "/" + std::to_string(ndf) + " = " + std::to_string(chi2_ndf);
    return Form("#chi^{2}/NDF = %.2f/%.2f = %.2f", chi2, ndf, chi2_ndf);
}

void plot_histogram::draw_th1d_result(TH1D* hist) {
    LOG_DEBUG.source("plot_histogram::draw_th1d_result") << "Start print th1d result";
    LOG_INFO.source("plot_histogram::drwa_th1d_result") << "Hist name: " << hist->GetName();

    TLatex* latex = new TLatex();
    set_tlatexStyle(latex);
    TLegend* legend = new TLegend(0.15, 0.7, 0.22, 0.8);

    std::string entry_text, mean_text, rms_text;
    double entry, mean, rms;
    double mean_error, rms_error;

    entry = hist->GetEntries();
    mean = hist->GetMean();
    rms = hist->GetRMS();
    mean_error = hist->GetMeanError();
    rms_error = hist->GetRMSError();

    legend->SetTextSize(0.04);
    legend->AddEntry(hist, "Data", "p");
    legend->Draw();

    entry_text = Form("Entry = %.0f", entry);
    mean_text = Form("Mean = %.2f #pm %.2f", mean, mean_error);
    rms_text = Form("RMS = %.2f #pm %.2f", rms, rms_error);
    latex->DrawLatexNDC(0.63, 0.84, entry_text.c_str());
    latex->DrawLatexNDC(0.63, 0.78, mean_text.c_str());
    latex->DrawLatexNDC(0.63, 0.72, rms_text.c_str());
}

void plot_histogram::print_fitResult(TF1* fit) {
    for(int i=0; i<fit->GetNpar();i++) {
        LOG_INFO.source("plot_histogram::print_TFitResult") << fit->GetParName(i) << ": " 
                                                            << fit->GetParameter(i) << " +/- " << fit->GetParError(i);
    }
}

void plot_histogram::draw_gausFit_result(TH1D* hist, TF1* fit) {
    LOG_DEBUG.source("plot_histogram::draw_gausFit_result") << "Start draw_gausFit_result";
    LOG_INFO.source("plot_histogram::draw_gausFit_result") << "Hist name: " << hist->GetName();
    LOG_INFO.source("plot_histogram::draw_gausFit_result") << "Fit name: " << fit->GetName();
    
    TLatex* latex = new TLatex();
    set_tlatexStyle(latex);
    TLegend* legend = new TLegend(0.15, 0.7, 0.22, 0.8);

    std::string mean_text, sigma_text, chi2ndf_text;
    double entry, mean, sigma, chi2, ndf;
    double mean_error, sigma_error;

    entry = hist->GetEntries();
    mean = fit->GetParameter(1);
    sigma = fit->GetParameter(2);
    mean_error = fit->GetParError(1);
    sigma_error = fit->GetParError(2);
    chi2 = fit->GetChisquare();
    ndf = fit->GetNDF();

    legend->SetTextSize(0.04);
    legend->AddEntry(hist, "Data", "p");
    legend->AddEntry(fit, "Fit", "l");
    legend->Draw();

    mean_text = get_tf1Parameter(fit, 1);
    sigma_text = get_tf1Parameter(fit, 2);
    chi2ndf_text = get_chi2ndfParameter(fit);
    latex->DrawLatexNDC(0.63, 0.84, Form("Entry = %.0f", entry));
    latex->DrawLatexNDC(0.63, 0.78, mean_text.c_str());
    latex->DrawLatexNDC(0.63, 0.72, sigma_text.c_str());
    latex->DrawLatexNDC(0.63, 0.66, chi2ndf_text.c_str());
}

void plot_histogram::draw_langausFit_result(TH1D* hist, TF1* fit) {
    LOG_DEBUG.source("plot_histogram::draw_langauFit_result") << "Start print_langausFit_result";
    LOG_INFO.source("plot_histogram::draw_langauFit_result") << "Hist name: " << hist->GetName();
    LOG_INFO.source("plot_histogram::draw_langauFit_result") << "Fit name: " << fit->GetName();
    
    TLatex* latex = new TLatex();
    set_tlatexStyle(latex);
    TLegend* legend = new TLegend(0.15, 0.7, 0.22, 0.8);

    std::string width_text, mp_text, area_text, gsigma_text, chi2ndf_text;
    double entry, width, mp, area, gsigma, chi2, ndf;
    double width_error, mp_error, area_error, gsigma_error;

    entry = hist->GetEntries();
    width = fit->GetParameter(0);
    mp = fit->GetParameter(1);
    area = fit->GetParameter(2);
    gsigma = fit->GetParameter(3);
    width_error = fit->GetParError(0);
    mp_error = fit->GetParError(1);
    area_error = fit->GetParError(2);
    gsigma_error = fit->GetParError(3);
    chi2 = fit->GetChisquare();
    ndf = fit->GetNDF();

    legend->SetTextSize(0.04);
    legend->AddEntry(hist, "Data", "p");
    legend->AddEntry(fit, "Fit", "l");
    legend->Draw();

    width_text = get_tf1Parameter(fit, 0);
    mp_text = get_tf1Parameter(fit, 1);
    area_text = get_tf1Parameter(fit, 2);
    gsigma_text = get_tf1Parameter(fit, 3);
    chi2ndf_text = get_chi2ndfParameter(fit);
    latex->DrawLatexNDC(0.63, 0.84, Form("Entry = %.0f", entry));
    latex->DrawLatexNDC(0.63, 0.78, width_text.c_str());
    latex->DrawLatexNDC(0.63, 0.72, mp_text.c_str());
    latex->DrawLatexNDC(0.63, 0.66, area_text.c_str());
    latex->DrawLatexNDC(0.63, 0.60, gsigma_text.c_str());
    latex->DrawLatexNDC(0.63, 0.54, chi2ndf_text.c_str());
}

// estimate FWHM of histogram for TH1D
std::pair<double, double> plot_histogram::estimate_fwhm(TH1D* hist) {
    double mean, peak, rms;
    double half_left, half_right;
    double center, fwhm;
    
    mean = hist->GetMean();
    peak = hist->GetMaximum();
    rms = hist->GetRMS();
    half_left = hist->FindFirstBinAbove(peak/2);
    half_right = hist->FindLastBinAbove(peak/2);
    center = (hist->GetBinCenter(half_left) + hist->GetBinCenter(half_right))/2;
    fwhm = hist->GetBinCenter(half_right) - hist->GetBinCenter(half_left);
    if (fwhm < 2 * hist->GetBinWidth(1)) {
        //print_message("FWHM is less than 2 times bin width", RED);
        LOG_WARNING.source("plot_histogram::estimate_fwhm()") << "FWHM is less than 2 times bin width.";
        return {rms, mean};
    }
    return {fwhm, center};
}

// optimise fitting function with Landau + Gaussian from langaus.C and set TF1 style
TF1* plot_histogram::optimise_hist_langaus(TH1D* hist, int color) {
    double area;
    std::vector<double> fit_range;
    std::vector<double> fwhm_pars, center_pars, area_pars;
    std::string fitname;
    std::vector<double> startvals, parlimitlow, parlimithigh;

    auto [fwhm, center] = estimate_fwhm(hist);
    fit_range = {0.3 * hist->GetMean(), 1.5 * hist->GetMean()};
    // fit_range_min = 0.3 * hist->GetMean();
    // fit_range_max = 1.5 * hist->GetMean();
    area = hist->Integral();

    fwhm_pars = {fwhm * 0.1, fwhm * 0, fwhm * 1};
    center_pars ={center, center * 0.5, center * 2};
    area_pars = {area * 10, area * 5, area * 100};

    fitname = std::string("fLangaus_") + hist->GetName();
    TF1* fit = new TF1(fitname.c_str(), langaufun, fit_range[0], fit_range[1], 4);
    set_tf1Style(fit, color);
    startvals = {fwhm_pars[0], center_pars[0], area_pars[0], fwhm_pars[0]};
    parlimitlow = {fwhm_pars[1], center_pars[1], area_pars[1], fwhm_pars[1]};
    parlimithigh = {fwhm_pars[2], center_pars[2], area_pars[2], fwhm_pars[2]};
    fit->SetParameters(startvals.data());
    fit->SetParNames("Width", "MP", "Area", "GSigma"); // MP is Most Probable Value, Area is integral of Landau, GSigma is gaussian sigma
    for (int i = 0; i < 4; ++i) {
        fit->SetParLimits(i, parlimitlow[i], parlimithigh[i]);
    }

    LOG_INFO.source("plot_histogram::optimise_hist_langaus") << "Hist name is " << hist->GetName() << "/Fit name is " << fitname;
    LOG_INFO.source("plot_histogram::optimise_hist_langaus") << "Fit range: " << fit_range[0] << " -> " << fit_range[1]; 
    hist->Fit(fit, "RLQ", "", fit_range[0], fit_range[1]);
    print_fitResult(fit);
    LOG_INFO.source("plot_histogram::optimise_hist_langaus") << "Chi2/ndf: " << fit->GetChisquare() << "/" << fit->GetNDF();
    if(fit->GetChisquare() / fit->GetNDF() > 10) {
        //print_message("Chi2/ndf is large, check the fitting!", RED);
        LOG_WARNING.source("plot_histogram::optimise_hist_langaus()") << "Chi2/ndf in [" << hist->GetName() << "] is larger than 10.";
    }
    return fit;
}

// optimise fitting function with Gaussian and set TF1 style
TF1* plot_histogram::optimise_hist_gaus(TH1D* hist, int color) {
    double peak, mean, rms;
    double half_left, half_right;
    double center, fwhm;
    double fit_range, fit_range_min, fit_range_max;
    std::string fit_name;

    peak = hist->GetMaximum();
    mean = hist->GetMean();
    rms = hist->GetRMS();
    half_left = hist->FindFirstBinAbove(peak/2);
    half_right = hist->FindLastBinAbove(peak/2);
    center = (hist->GetBinCenter(half_left) + hist->GetBinCenter(half_right))/2;
    fwhm = hist->GetBinCenter(half_right) - hist->GetBinCenter(half_left);

    fit_name = std::string("fGaus_") + hist->GetName();
    fit_range = std::min(2 * rms, 2 * fwhm);
    fit_range_min = center - fit_range;
    fit_range_max = center + fit_range;

    LOG_INFO.source("plot_histogram::optimise_hist_gaus") << "Hist name is " << hist->GetName() << "/Fit name is " << fit_name;
    LOG_INFO.source("plot_histogram::optimise_hist_gaus") << "Fit range: " << fit_range_min << " -> " << fit_range_max; 
    
    TF1* fit = new TF1(fit_name.c_str(), "gaus", -60, 60);
    set_tf1Style(fit, color);
    hist->Fit(fit, "RLQ", "", fit_range_min, fit_range_max);
    print_fitResult(fit);
    LOG_INFO.source("plot_histogram::optimise_hist_gaus") << "Chi2/ndf: " << fit->GetChisquare() << "/" << fit->GetNDF();
    if(fit->GetChisquare() / fit->GetNDF() > 10) {
        LOG_WARNING.source("plot_histogram::optimise_hist_gaus()") << "Chi2/ndf in [" << hist->GetName() << "] is larger than 10.";
    }
    return fit;
}

TF1* plot_histogram::optimise_hist_doublegaus(TH1D* hist, int color) {
    double peak, mean, rms;
    double x1, x2, x3, x4;
    double mean1, mean2, sigma1, sigma2;
    double fit_range_min, fit_range_max;
    std::string fit_name;

    peak = hist->GetMaximum();
    mean = hist->GetMean();
    rms = hist->GetRMS();

    x1 = hist->GetBinCenter(hist->FindFirstBinAbove(peak / 2));  // 左端
    x2 = x1 + rms/4;  // 左2点目;  // 右2点目
    x4 = hist->GetBinCenter(hist->FindLastBinAbove(peak / 2));  // 右端
    x3 = x4 - rms/4;  // 右2点目

    mean1 = (x1 + x2) / 2.0;
    mean2 = (x3 + x4) / 2.0;

    sigma1 = (x2 - x1) / 2.355;
    sigma2 = (x4 - x3) / 2.355;

    fit_range_min = mean1 - 3 * sigma1;
    fit_range_max = mean2 + 3 * sigma2;

    fit_name = std::string("fDoubleGaus_") + hist->GetName();

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << BLUE << "Hist name: " << hist->GetName() << RESET << std::endl;
    std::cout << BLUE << "Fit name: " << fit_name << RESET << std::endl;
    std::cout << "Fit range: " << fit_range_min << " -> " << fit_range_max << std::endl;

    TF1* fit = new TF1(fit_name.c_str(), twogausfun, 0, 10000, 6);

    fit->SetParNames("A1", "Mean1", "Sigma1", "A2", "Mean2", "Sigma2");

    //fit->SetParameters(peak, mean1, sigma1, peak / 2, mean2, sigma2);
    fit->SetParameters(190, -3, 0.8, 155, 4, 0.8);

    set_tf1Style(fit, color);


    int fitStatus = hist->Fit(fit, "RL", "", fit_range_min, fit_range_max);

    if (fitStatus != 0) {
        std::cerr << "Error: Fit failed with status " << fitStatus << std::endl;
        return nullptr;
    }

    std::cout << "Chi2/ndf: " << fit->GetChisquare() << "/" << fit->GetNDF() << std::endl;
    if (fit->GetChisquare() / fit->GetNDF() > 10) {
        print_message("Chi2/ndf is large, check the fitting!", RED);
    }
    std::cout << "-------------------------------------------------------------------------" << std::endl;

    return fit;
}

// draw histogram from clusteringSpatial module as reference plots
void plot_histogram::draw_ref_clusteringSpatial(const std::string& refname, TFile* input_file, TFile* output_file) {
    LOG_DEBUG.source("plot_histogram::draw_ref_clusteringSpatial") << "Start draw_ref_clusteringSpatial";
    LOG_INFO.source("plot_histogram::draw_ref_clusteringSpatial") << "refname: " + refname;
    
    input_dirName = "ClusteringSpatial/" + refname + "/";
    output_dirName = "Hitmap/" + refname + "/";

    LOG_STATUS.source("plot_histogram::draw_ref_clusteringSpatial") << "Inport histogram from " << input_dirName;
    hHitmapLocal = (TH2D*)input_file->Get((input_dirName + "clusterPositionLocal").c_str());
    hHitmapGlobal = (TH2D*)input_file->Get((input_dirName + "clusterPositionGlobal").c_str());
    if (hHitmapLocal == nullptr || hHitmapGlobal == nullptr) {
        LOG_ERROR.source("plot_histogram::draw_ref_clusteringSpatial") << "Histogram is not found at " << input_dirName;
        throw std::invalid_argument("Histogram is not found");
    }
    hProjectionX = hHitmapLocal->ProjectionX("hProjectionX");
    hProjectionY = hHitmapLocal->ProjectionY("hProjectionY");

    set_th2dStyle(hHitmapLocal);
    set_th2dStyle(hHitmapGlobal);
    set_th1dStyle(hProjectionX, kRed);
    set_th1dStyle(hProjectionY, kRed);

    LOG_STATUS.source("plot_histogram::draw_ref_clusteringSpatial") << "Optimise fitting";
    fProjectionX = optimise_hist_gaus(hProjectionX, kBlue);
    fProjectionY = optimise_hist_gaus(hProjectionY, kBlue);

    LOG_STATUS.source("plot_histogram::draw_ref_clusteringSpatial") << "Write histogram to " << output_dirName;
    output_file->cd(output_dirName.c_str());
    write_plot(hHitmapLocal, "hHitmapLocal");
    write_plot(hHitmapGlobal, "hHitmapGlobal");
    write_plot(hProjectionX, "hProjectionX");
    write_plot(hProjectionY, "hProjectionY");
    output_file->cd();
}

// draw histogram from clusteringAnalog module as DUT plots
void plot_histogram::draw_dut_clusteringAnalog(const std::string& dutname, TFile* input_file, TFile* output_file) {
    LOG_DEBUG.source("plot_histogram::draw_dut_clusteringAnalog") << "Start draw_dut_clusteringAnalog";
    LOG_INFO.source("plot_histogram::draw_dut_clusteringAnalog") << "dutname: " + dutname;

    input_dirName = "ClusteringAnalog/" + dutname + "/";
    output_dirName = "Hitmap/" + dutname + "/";

    LOG_STATUS.source("plot_histogram::draw_dut_clusteringAnalog") << "Inport histogram from " << input_dirName;
    hHitmapLocal = (TH2D*)input_file->Get((input_dirName + "clusterPositionLocal").c_str());
    hHitmapGlobal = (TH2D*)input_file->Get((input_dirName + "clusterPositionGlobal").c_str());
    if (hHitmapLocal == nullptr || hHitmapGlobal == nullptr) {
        LOG_ERROR.source("plot_histogram::draw_dut_clusteringAnalog") << "Histogram is not found at " << input_dirName;
        throw std::invalid_argument("Histogram is not found");
    }
    hProjectionX = hHitmapLocal->ProjectionX("hProjectionX");
    hProjectionY = hHitmapLocal->ProjectionY("hProjectionY");

    set_th2dStyle(hHitmapLocal);
    set_th2dStyle(hHitmapGlobal);
    set_th1dStyle(hProjectionX, kRed);
    set_th1dStyle(hProjectionY, kRed);

    LOG_STATUS.source("plot_histogram::draw_dut_clusteringAnalog") << "Optimise fitting";
    fProjectionX = optimise_hist_gaus(hProjectionX, kBlue);
    fProjectionY = optimise_hist_gaus(hProjectionY, kBlue);

    LOG_STATUS.source("plot_histogram::draw_dut_clusteringAnalog") << "Write histogram to " << output_dirName;
    output_file->cd(output_dirName.c_str());
    write_plot(hHitmapLocal, "hHitmapLocal");
    write_plot(hHitmapGlobal, "hHitmapGlobal");
    write_plot(hProjectionX, "hProjectionX");
    write_plot(hProjectionY, "hProjectionY");
    output_file->cd();
}

// draw correlation histogram from Correlations module
void plot_histogram::draw_correlation(const std::string& detector, TFile* input_file, TFile* output_file) {
    std::cout << "=========================================================================" << std::endl;
    print_message("Start draw_correlation", GREEN);
    print_message("detector: " + detector, GREEN);
    
    input_dirName = "Correlations/" + detector + "/";
    output_dirName = "Correlation/" + detector + "/";

    std::cout << "Inport histogram from " << input_dirName << std::endl;
    hCorrelationX = (TH1D*)input_file->Get((input_dirName + "correlationX").c_str());
    hCorrelationY = (TH1D*)input_file->Get((input_dirName + "correlationY").c_str());
    hCorrelationX_2D = (TH2D*)input_file->Get((input_dirName + "correlationX_2Dlocal").c_str());
    hCorrelationY_2D = (TH2D*)input_file->Get((input_dirName + "correlationY_2Dlocal").c_str());
    if (hCorrelationX == nullptr || hCorrelationY == nullptr || hCorrelationX_2D == nullptr || hCorrelationY_2D == nullptr) {
        print_message(Form("Histogram is not found at %s", input_dirName.c_str()), RED);
        return;
    }

    set_th1dStyle(hCorrelationX, kRed);
    set_th1dStyle(hCorrelationY, kRed);
    set_th2dStyle(hCorrelationX_2D);
    set_th2dStyle(hCorrelationY_2D);

    std::cout << "Optimise fitting" << std::endl;
    fCorrelationX = optimise_hist_gaus(hCorrelationX, kBlue);
    fCorrelationY = optimise_hist_gaus(hCorrelationY, kBlue);

    std::cout << "Write histogram to " << output_dirName << std::endl;
    output_file->cd(output_dirName.c_str());
    write_plot(hCorrelationX, "correlationX");
    write_plot(hCorrelationY, "correlationY");
    write_plot(hCorrelationX_2D, "correlationX_2D");
    write_plot(hCorrelationY_2D, "correlationY_2D");
    output_file->cd();

    std::cout << "Finish draw_correlation" << std::endl;
    std::cout << "detector: " << detector << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
}

// draw histogram from Tracking4D module
void plot_histogram::draw_Tracking4D(TFile* input_file, TFile* output_file) {
    std::cout << "=========================================================================" << std::endl;
    print_message("Draw Tracking4D", GREEN);

    std::string input_dirName, output_dirName;

    input_dirName = "Tracking4D/";
    output_dirName = "Track";

    std::cout << "Inport histogram from " << input_dirName << std::endl;
    hTrackChi2 = (TH1D*)input_file->Get((input_dirName + "trackChi2").c_str());
    hTrackChi2ndof = (TH1D*)input_file->Get((input_dirName + "trackChi2ndof").c_str());
    hTrackAngleX = (TH1D*)input_file->Get((input_dirName + "trackAngleX").c_str());
    hTrackAngleY = (TH1D*)input_file->Get((input_dirName + "trackAngleY").c_str());
    if (hTrackChi2 == nullptr || hTrackChi2ndof == nullptr || hTrackAngleX == nullptr || hTrackAngleY == nullptr) {
        print_message(Form("Histogram is not found at %s", input_dirName.c_str()), RED);
        return;
    }

    set_th1dStyle(hTrackChi2, kRed);
    set_th1dStyle(hTrackChi2ndof, kRed);
    set_th1dStyle(hTrackAngleX, kRed);
    set_th1dStyle(hTrackAngleY, kRed);

    std::cout << "Optimise fitting" << std::endl;
    fTrackAngleX = optimise_hist_gaus(hTrackAngleX, kBlue);
    fTrackAngleY = optimise_hist_gaus(hTrackAngleY, kBlue);

    std::cout << "Write histogram to " << output_dirName << std::endl;
    output_file->cd(output_dirName.c_str());
    write_plot(hTrackChi2, "trackChi2");
    write_plot(hTrackChi2ndof, "trackChi2ndof");
    write_plot(hTrackAngleX, "trackAngleX");
    write_plot(hTrackAngleY, "trackAngleY");
    output_file->cd();

    std::cout << "Finish draw_Tracking4D" << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
}

// draw histogram from some modules (ClusteringSpatial, Tracking4D) as reference plots
void plot_histogram::draw_ref_analysis(std::string refname, TFile* input_file, TFile* output_file) {
    std::cout << "=========================================================================" << std::endl;
    print_message("Start draw_ref_analysis", GREEN);
    print_message("refname: " + refname, GREEN);

    std::string input_dirName_ClSpatial, input_dirName_Tg4D, output_dirName;

    input_dirName_ClSpatial = "ClusteringSpatial/" + refname + "/";
    input_dirName_Tg4D = "Tracking4D/" + refname + "/local_residuals/";
    output_dirName = "AnalysisRef/" + refname + "/";

    std::cout << "Inport histogram from " << input_dirName_ClSpatial << std::endl;
    hClusterSize = (TH1D*)input_file->Get((input_dirName_ClSpatial + "clusterSize").c_str());
    hClusterSeedCharge = (TH1D*)input_file->Get((input_dirName_ClSpatial + "clusterSeedCharge").c_str());
    hClusterCharge = (TH1D*)input_file->Get((input_dirName_ClSpatial + "clusterCharge").c_str());
    hClusterMultiplicity = (TH1D*)input_file->Get((input_dirName_ClSpatial + "clusterMultiplicity").c_str());
    hResidualsX = (TH1D*)input_file->Get((input_dirName_Tg4D + "LocalResidualsX").c_str());
    hResidualsY = (TH1D*)input_file->Get((input_dirName_Tg4D + "LocalResidualsY").c_str());
    if (hClusterSize == nullptr || hClusterSeedCharge == nullptr || hClusterCharge == nullptr || hClusterMultiplicity == nullptr || hResidualsX == nullptr || hResidualsY == nullptr) {
        print_message(Form("Histogram is not found at %s", input_dirName_ClSpatial.c_str()), RED);
        return;
    }

    set_th1dStyle(hClusterSize, kRed);
    set_th1dStyle(hClusterSeedCharge, kRed);
    set_th1dStyle(hClusterCharge, kRed);
    set_th1dStyle(hClusterMultiplicity, kRed);
    set_th1dStyle(hResidualsX, kRed);
    set_th1dStyle(hResidualsY, kRed);

    std::cout << "Optimise fitting" << std::endl;
    fResidualsX = optimise_hist_gaus(hResidualsX, kBlue);
    fResidualsY = optimise_hist_gaus(hResidualsY, kBlue);

    std::cout << "Write histogram to " << output_dirName << std::endl;
    output_file->cd(output_dirName.c_str());
    write_plot(hClusterSize, "clusterSize");
    write_plot(hClusterSeedCharge, "clusterSeedCharge");
    write_plot(hClusterCharge, "clusterCharge");
    write_plot(hClusterMultiplicity, "clusterMultiplicity");
    write_plot(hResidualsX, "residualsX");
    write_plot(hResidualsY, "residualsY");
    output_file->cd();

    std::cout << "Finish draw_ref_analysis" << std::endl;
    std::cout << "refname: " << refname << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
}

// draw histogram from AnalysisCE65 module as DUT plots
void plot_histogram::draw_AnalysisCE65(std::string dutname, TFile* input_file, TFile* output_file) {
    std::cout << "=========================================================================" << std::endl;
    print_message("Start draw_AnalysisCE65", GREEN);
    print_message("dutname: " + dutname, GREEN);

    input_dirName = "AnalysisCE65/" + dutname + "/";
    output_dirName = "AnalysisCE65/" + dutname + "/";

    std::cout << "Inport histogram from " << input_dirName << std::endl;
    hCutHisto = (TH1D*)input_file->Get((input_dirName + "hCutHisto").c_str());
    hClusterSize = (TH1D*)input_file->Get((input_dirName + "cluster/clusterSize").c_str());
    hClusterSize_SeedCut = (TH1D*)input_file->Get((input_dirName + "cluster/clusterShape_SeedCut").c_str());
    hClusterSeedCharge = (TH1D*)input_file->Get((input_dirName + "cluster/clusterSeedCharge").c_str());
    hClusterCharge = (TH1D*)input_file->Get((input_dirName + "cluster/clusterCharge").c_str());
    hClusterNeighborCharge = (TH1D*)input_file->Get((input_dirName + "cluster/clusterNeighborsCharge").c_str());
    hClusterNeighborChargeSum = (TH1D*)input_file->Get((input_dirName + "cluster/clusterNeighborsChargeSum").c_str());
    hClusterPositionLocal = (TH2D*)input_file->Get((input_dirName + "cluster/clusterPositionLocal").c_str());
    hResidualsX = (TH1D*)input_file->Get((input_dirName + "local_residuals/residualsX").c_str());
    hResidualsY = (TH1D*)input_file->Get((input_dirName + "local_residuals/residualsY").c_str());
    if (hCutHisto == nullptr || hClusterSize == nullptr || hClusterSize_SeedCut == nullptr || hClusterSeedCharge == nullptr ||
        hClusterCharge == nullptr || hClusterNeighborCharge == nullptr || hClusterNeighborChargeSum == nullptr ||
        hClusterPositionLocal == nullptr || hResidualsX == nullptr || hResidualsY == nullptr) {
        print_message(Form("Histogram is not found at %s", input_dirName.c_str()), RED);
        return;
    }

    set_th1dStyle(hCutHisto, kRed);
    set_th1dStyle(hClusterSize, kRed);
    set_th1dStyle(hClusterSize_SeedCut, kRed);
    set_th1dStyle(hClusterSeedCharge, kRed);
    set_th1dStyle(hClusterCharge, kRed);
    set_th1dStyle(hClusterNeighborCharge, kRed);
    set_th1dStyle(hClusterNeighborChargeSum, kRed);
    set_th2dStyle(hClusterPositionLocal);
    set_th1dStyle(hResidualsX, kRed);
    set_th1dStyle(hResidualsY, kRed);

    std::cout << "Optimise fitting" << std::endl;
    fClusterCharge = optimise_hist_langaus(hClusterCharge, kBlue);
    fResidualsX = optimise_hist_gaus(hResidualsX, kBlue);
    fResidualsY = optimise_hist_gaus(hResidualsY, kBlue);

    std::cout << "Write histogram to " << output_dirName << std::endl;
    output_file->cd(output_dirName.c_str());
    write_plot(hCutHisto, "CutHisto");
    write_plot(hClusterSize, "clusterSize");
    write_plot(hClusterSize_SeedCut, "clusterSize_SeedCut");
    write_plot(hClusterSeedCharge, "clusterSeedCharge");
    write_plot(hClusterCharge, "clusterCharge");
    write_plot(hClusterNeighborCharge, "clusterNeighborCharge");
    write_plot(hClusterNeighborChargeSum, "clusterNeighborChargeSum");
    write_plot(hClusterPositionLocal, "clusterPositionLocal");
    write_plot(hResidualsX, "residualsX");
    write_plot(hResidualsY, "residualsY");
    output_file->cd();

    std::cout << "Finish draw_AnalysisCE65" << std::endl;
    std::cout << "dutname: " << dutname << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
}

// run all functions
void plot_histogram::run() {
    set_rootStyle();

    refnames = {"ALPIDE_0", "ALPIDE_1", "ALPIDE_2", "ALPIDE_4", "ALPIDE_5", "ALPIDE_6"};
    dutnames = {"CE65_3"};
    input_fileDir = "/home/towa/alice3/hist/kek_sim/";
    base_filename = "05_analysis_kek_sim_CE65_p15_std_6ALPIDE_202502131448.root";
    input_filename = input_fileDir + base_filename;
    output_filename = "plot_" + base_filename;
    rootFile_dir_names = {"Hitmap", "Correlation", "AnalysisRef", "Track", "AnalysisCE65"};

    input_file = new TFile(input_filename.c_str());
    output_file = new TFile(output_filename.c_str(), "RECREATE");

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
}

// // constructor
// plot_histogram::plot_histogram() {
//     std::cout << "plot_histogram object is created" << std::endl;
// }

void plot_histogram::set_gausHist_range(TH1D* hist, int min, int max) {
    hist->GetXaxis()->SetRangeUser(min, max);
}

void plot_histogram::set_langausHist_range(TH1D* hist, int min, int max) {
    hist->GetXaxis()->SetRangeUser(min, max);
}

void plot_histogram::search_gausHist_range(TH1D* hist, double threshold) {
    int bin_min, bin_max;
    int bin;

    bin_min = hist->FindFirstBinAbove(threshold);
    bin_max = hist->FindLastBinAbove(threshold);
    bin = (bin_max - bin_min) / 2; // round down

    hist->GetXaxis()->SetRangeUser(-bin, bin);
}

void plot_histogram::search_langausHist_range(TH1D* hist, double threshold) {
    int bin_min, bin_max;

    //bin_min = hist->FindFirstBinAbove(threshold);
    bin_min = 0;
    bin_max = hist->FindLastBinAbove(threshold);

    hist->GetXaxis()->SetRangeUser(bin_min, bin_max + 1);
}

void plot_histogram::set_TH1D_position(TCanvas* c){
    c->SetTopMargin(0.04);
    c->SetBottomMargin(0.16);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.03);
}

void plot_histogram::set_TH1D_label(TH1D* hist_1D, double x_range, double y_range){
    hist_1D->SetStats(0);
    hist_1D->GetXaxis()->SetRangeUser(1, x_range);
    hist_1D->GetYaxis()->SetRangeUser(0, y_range);
    hist_1D->GetXaxis()->SetTitleOffset(0.8);
    hist_1D->GetYaxis()->SetTitleOffset(0.7);
    hist_1D->SetTitleSize(0.08, "XY");
    hist_1D->SetLabelSize(0.06, "XY");
}

void plot_histogram::set_TH1D_point(TH1D* hist_1D, int color, int style){
    hist_1D->SetMarkerColor(color);
    hist_1D->SetLineColor(color);
    hist_1D->SetMarkerStyle(style);
    hist_1D->SetMarkerSize(1.2);
}

void plot_histogram::makeDirectory(TFile *file, const std::vector<std::string> &dirnames) {
    for (int i = 0; i < dirnames.size(); ++i) {
        TDirectory *dir = file->mkdir(dirnames[i].c_str());
    }
    file->cd();
}

void plot_histogram::moveToDirectory(TFile *file, const std::string &dirName) {
    TDirectory *dir = (TDirectory *)file->Get(dirName.c_str());
    if (dir) {
        dir->cd();
    } else {
        std::cerr << "Error: Directory " << dirName << " does not exist!" << std::endl;
    }
}

void plot_histogram::write_canvas_Ref(TH1D* hist_x, TF1* fit_x, TH1D* hist_y, TF1* fit_y) {
    TCanvas *c_x = new TCanvas("c", "c", 800, 600);
    set_TH1D_position(c_x);
    hist_x->Draw("PE");
    fit_x->Draw("sameL");
    c_x->Write("ProjectionX_fit");
    TCanvas *c_y = new TCanvas("c", "c", 800, 600);
    set_TH1D_position(c_y);
    hist_y->Draw("PE");
    fit_y->Draw("sameL");
    c_y->Write("ProjectionY_fit");
}

void plot_histogram::saveCanvasesToPDF(const char* rootFileName, const char* dirName, const char* outputPDF) {
    TFile *file = TFile::Open(rootFileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file " << rootFileName << std::endl;
        return;
    }

    TDirectory *dir = dynamic_cast<TDirectory*>(file->Get(dirName));
    if (!dir) {
        std::cerr << "Error: Directory " << dirName << " not found in ROOT file." << std::endl;
        file->Close();
        return;
    }

    std::string pdfFile = outputPDF;
    std::string pdfStart = pdfFile + "[";
    std::string pdfEnd = pdfFile + "]";

    bool firstCanvas = true;

    TIter next(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        if (std::string(key->GetClassName()) == "TCanvas") {
            TCanvas *canvas = dynamic_cast<TCanvas*>(key->ReadObj());
            if (!canvas) continue;

            if (firstCanvas) {
                canvas->Print(pdfStart.c_str());
                firstCanvas = false;
            }

            canvas->Print(pdfFile.c_str());
        }
    }

    if (!firstCanvas) {
        TCanvas dummy;
        dummy.Print(pdfEnd.c_str());
    } else {
        std::cerr << "No TCanvas objects found in directory " << dirName << std::endl;
    }

    file->Close();
    delete file;
}

/// or filtering
std::vector<std::string> plot_histogram::or_filter_filenames(std::vector<std::string> filenames, std::vector<std::string> filter) {
    std::vector<std::string> filtered_filenames;
    for(int i=0; i<filenames.size(); ++i) {
        for(int j=0; j<filter.size(); ++j) {
            if(filenames[i].find(filter[j]) != std::string::npos) {
                filtered_filenames.push_back(filenames[i]);
            }
        }
    }
    //std::sort(filtered_filenames.begin(), filtered_filenames.end());
    return filtered_filenames;
}

// and filtering
std::vector<std::string> plot_histogram::and_filter_filenames(std::vector<std::string> filenames, std::vector<std::string> filter) {
    
    std::vector<std::string> filtered_filenames;
    
    for (const auto& filename : filenames) {
        bool all_found = true;
        
        for (const auto& f : filter) {
            if (filename.find(f) == std::string::npos) {  
                all_found = false;
                break;
            }
        }

        if (all_found) {
            filtered_filenames.push_back(filename);
        }
    }
    //std::sort(filtered_filenames.begin(), filtered_filenames.end());
    return filtered_filenames;
}