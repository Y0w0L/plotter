#ifndef PLOT_HISTOGRAM_HPP
#define PLOT_HISTOGRAM_HPP

#include <TFile.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TError.h>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <variant>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <TKey.h>

#include "Messenger/Messenger.h"
//#include "langaus.C"

class plot_histogram {
public:
    /**
    * @brief Constructor for this plot_histogram object
    **/
    plot_histogram();
    ~plot_histogram();
    // function
    /**
    * @brief Set ROOT style for the plot_histogram object
    **/
    virtual void set_rootStyle();
    /**
    * @brief Set TH1D histogram style for the plot_histogram object
    * @param hist Pointer to TH1D object
    * @param color Marker color
    **/
    virtual void set_th1dStyle(TH1D* hist, int color);
    /**
    * @brief Set TH2D histogram style for the plot_histogram object
    * @param hist Pointer to TH2D object
    **/
    virtual void set_th2dStyle(TH2D* hist);
    /**
    * @brief Set TF1 style for the plot_histogram object
    * @param fit Pointer to TF1 object
    * @param color Marker color
    **/
    static void set_tf1Style(TF1* fit, int color);
    /**
    * @brief Set TLatex style for the plot_histogram object
    * @param latex Pointer to TLatex object
    **/
    static void set_tlatexStyle(TLatex* latex);
    /**
    * @brief Remove white space in the string
    * @param word String to remove white space
    **/
    static void remove_whiteSpace(std::string& word);
    /**
    * @brief Replace string in the string
    * @param word String to replace
    * @param target Target word
    * @param replacement Replacement word
    **/
    static void replace_string(std::string& word,
                               const std::string& target,
                               const std::string& replacement);
    /**
    * @brief Search word in the string
    * @param word String to search
    * @param target Target word
    **/
    static bool search_word(const std::string& word, const std::string& target);
    /**
    * @brief Print message with color
    * @param message Message to print
    * @param color Color of the message
    **/
    static void print_message(std::string message, const char* color);
    /**
    * @brief Format uncertainty value
    * @param value Value
    * @param error Error
    **/
    static std::string format_uncertainty(double value, double error);
    /**
    * @brief Write plot to canvas
    * @param hist Pointer to TH1D or TH2D object
    * @param hist_name Name of the histogram
    **/
    virtual void write_plot(std::variant<TH1D*, TH2D*> hist, const std::string& hist_name);
    /**
    * @brief Get parameter of TF1 object
    * @param fit Pointer to TF1 object
    * @param par Parameter number
    **/
    static std::string get_tf1Parameter(TF1* fit, int par);
    /**
    * @brief Get chi2/ndf parameter of TF1 object
    * @param fit Pointer to TF1 object
    **/
    static std::string get_chi2ndfParameter(TF1* fit);
    /**
    * @brief Draw TH1D result
    * @param hist Pointer to TH1D object
    **/
    static void draw_th1d_result(TH1D* hist);
    /**
    * @brief Draw gaus fitting result
    * @param hist Pointer to TH1D object
    * @param fit Pointer to TF1 object
    **/
    static void draw_gausFit_result(TH1D* hist, TF1* fit);
    /**
    * @brief Draw landau gaussian fitting result // from langaus.C
    * @param hist Pointer to TH1D object
    * @param fit Pointer to TF1 object
    **/
    static void draw_langausFit_result(TH1D* hist, TF1* fit);
    /**
    * @brief Estimate FWHM of the histogram for TH1D
    * @param hist Pointer to TH1D object
    **/
    static std::pair<double, double> estimate_fwhm(TH1D* hist);
    /**
    * @brief Optimise fitting function with Landau-Gaussian
    * @param hist Pointer to TH1D object
    * @param color TF1 line color
    **/
    static TF1* optimise_hist_langaus(TH1D* hist, int color);
    /**
    * @brief Optimise fitting function with Gaussian
    * @param hist Pointer to TH1D object
    * @param color TF1 line color
    **/
    static TF1* optimise_hist_gaus(TH1D* hist, int color);
    TF1* optimise_hist_doublegaus(TH1D* hist, int color);
    /**
    * @brief Draw histogram from clusteringSpatial module as reference plots
    * @param refname Name of the reference module
    * @param input_file Pointer to the input ROOT file
    * @param output_file Pointer to the output ROOT file
    **/
    void draw_ref_clusteringSpatial(std::string refname, TFile* input_file, TFile* output_file);
    /**
    * @brief Draw histogram from clusteringAnalog module as DUT plots
    * @param dutname Name of the DUT module
    * @param input_file Pointer to the input ROOT file
    * @param output_file Pointer to the output ROOT file
    **/
    void draw_dut_clusteringAnalog(std::string dutname, TFile* input_file, TFile* output_file);
    /**
    * @brief Draw correlation histogram from Correlations module
    * @param detector Name of the detector module
    * @param input_file Pointer to the input ROOT file
    * @param output_file Pointer to the output ROOT file
    **/
    void draw_correlation(std::string detector, TFile* input_file, TFile* output_file);
    /**
    * @brief Draw histogram from Tracking4D module
    * @param input_file Pointer to the input ROOT file
    * @param output_file Pointer to the output ROOT file
    **/
    void draw_Tracking4D(TFile* input_file, TFile* output_file);
    /**
    * @brief Draw histogram from Correlations and ClusteringAnalog module as reference plots
    * @param refname Name of the reference module
    * @param input_file Pointer to the input ROOT file
    * @param output_file Pointer to the output ROOT file
    **/
    void draw_ref_analysis(std::string refname, TFile* input_file, TFile* output_file);
    /**
    * @brief Draw histogram from AnalysisCE65 module as DUT plots
    * @param dutname Name of the DUT module
    * @param input_file Pointer to the input ROOT file
    * @param output_file Pointer to the output ROOT file
    **/
    void draw_AnalysisCE65(std::string dutname, TFile* input_file, TFile* output_file);
    void saveCanvasesToPDF(const char* rootFileName, const char* dirName, const char* outputPDF);
    std::vector<std::string> or_filter_filenames(std::vector<std::string> filenames, std::vector<std::string> filter);
    std::vector<std::string> and_filter_filenames(std::vector<std::string> filenames, std::vector<std::string> filter);
    
    /**
    * @brief Run the plot_histogram object
    **/
    virtual void run();
    // histogram
    TH2D *hHitmapLocal, *hHitmapGlobal, *hClusterPositionLocal;
    TH1D *hProjectionX, *hProjectionY;
    TH1D *hCorrelationX, *hCorrelationY;
    TH2D *hCorrelationX_2D, *hCorrelationY_2D;
    TH1D *hTrackChi2, *hTrackChi2ndof, *hTrackAngleX, *hTrackAngleY;
    TH1D *hClusterSize, *hClusterSize_SeedCut, *hClusterSeedCharge, *hClusterCharge, *hClusterNeighborCharge,
         *hClusterNeighborChargeSum;
    TH1D* hCutHisto;
    TH1D* hClusterMultiplicity;
    TH1D* hResidualsX, *hResidualsY;
    // fitting
    TF1 *fProjectionX, *fProjectionY;
    TF1 *fCorrelationX, *fCorrelationY;
    TF1 *fTrackAngleX, *fTrackAngleY;
    TF1 *fResidualsX, *fResidualsY;
    TF1 *fClusterCharge;
    // console color
    #define RESET   "\033[0m"
    #define RED     "\033[31m"
    #define GREEN   "\033[32m"
    #define YELLOW  "\033[33m"
    #define BLUE    "\033[34m"
private:
    std::string name_;
    // function
    /**
    * @brief Set Gaussian range
    * @param hist Pointer to TH1D object
    * @param min Minimum value of the range
    * @param max Maximum value of the range
    **/
    void set_gausHist_range(TH1D* hist, int min, int max);
    /**
    * @brief Set Landau-Gaussian range
    * @param hist Pointer to TH1D object
    * @param min Minimum value of the range
    * @param max Maximum value of the range
    **/
    void set_langausHist_range(TH1D* hist, int min, int max);
    /**
    * @brief Search Gaussian histogram range
    * @param hist Pointer to TH1D object
    * @param threshold Threshold value
    **/
    void search_gausHist_range(TH1D* hist, double threshold);
    /**
    * @brief Search Landau-Gaussian histogram range
    * @param hist Pointer to TH1D object
    * @param threshold Threshold value
    **/
    void search_langausHist_range(TH1D* hist, double threshold);
    /**
    * @brief Set TH1D position
    * @param c Pointer to TCanvas object
    **/
    void set_TH1D_position(TCanvas* c);
    /**
    * @brief Set TH1D label
    * @param hist Pointer to TH1D object
    * @param x_range X-axis range
    * @param y_range Y-axis range
    **/
    void set_TH1D_label(TH1D* hist, double x_range, double y_range);
    /**
    * @brief Set TH1D point
    * @param hist Pointer to TH1D object
    * @param color Marker color
    * @param style Marker style
    **/
    void set_TH1D_point(TH1D* hist, int color, int style);
    /**
    * @brief Make directory in the ROOT file
    * @param file Pointer to TFile object
    * @param dirnames Vector of directory names
    **/
    void makeDirectory(TFile* file, const std::vector<std::string>& dirnames);
    /**
    * @brief Move to directory in the ROOT file
    * @param file Pointer to TFile object
    * @param dirName Name of the directory
    **/
    void moveToDirectory(TFile* file, const std::string& dirName);
    /**
    * @brief Write histogram to ROOT file
    * @param hist_x Pointer to TH1D object
    * @param fit_x Pointer to TF1 object
    * @param hist_y Pointer to TH1D object
    * @param fit_y Pointer to TF1 object
    **/
    void write_canvas_Ref(TH1D* hist_x, TF1* fit_x, TH1D* hist_y, TF1* fit_y);
    std::string input_dirName, output_dirName;
    std::string input_dirName_ClSpatial, input_dirName_Tg4D;
    std::vector<std::string> refnames, dutnames;
    std::vector<std::string> rootFile_dir_names;
    std::vector<TDirectory*> rootFile_dirs;
    std::string input_fileDir, base_filename;
    std::string input_filename, output_filename;
    TFile *input_file, *output_file;
    TDirectory *hitmap, *correlation, *analysisRef, *track, *analysisCE65;
};

#endif // PLOT_HISTOGRAM_HPP