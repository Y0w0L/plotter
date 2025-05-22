#include "plot_histogram/plot_histogram.hpp"

class plot_simulation : public plot_histogram {
    public:
        /**
        * @brief Constructor for this plot_simulation object
        **/
        plot_simulation();

        /**
        * @brief Set ROOT style for the plot_simulation object
        **/
        void set_rootStyle() override;

        /**
        * @brief Set TH1D histogram style for the plot_simulation object
        * @param hist Pointer to TH1D object
        * @param color Marker color
        **/
        void set_th1dStyle(TH1D* hist, int color);

        /**
        * @brief Set TH2D histogram style for the plot_simulation object
        * @param hist Pointer to TH2D object
        **/
        void set_th2dStyle(TH2D* hist);

        static TH1D* search_TH1D_inCanvas(TCanvas* c);

        /**
        * @brief Search TF1 object in canvas or TH1 object
        * @param c Pointer to TCanvas object
        **/
        static TF1* search_TF1_inCanvas(TCanvas* c);

        /**
        * @brief Write plot to canvas
        * @param hist Pointer to TH1D or TH2D object
        * @param hist_name Name of the histogram
        **/
        void write_plot(std::variant<TH1D*, TH2D*> hist, const std::string& hist_name);
        // static void tlatex_style(TLatex* latex);
        // static void draw_langausFit_result(TH1D* hist, TF1* fit);
        // static void draw_gausFit_result(TH1D* hist, TF1* fit);

        /**
        * @brief Run the plot_simulation object
        **/
        void run() override;

    private:
        /**
        * @brief Set cluster size distribution style
        * @param hist Pointer to TH1D object
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_clSizeStyle(TH1D* hist, double range, int rebin);

        /**
        * @brief Set Correlation distribution style
        * @param hist Pointer to TH1D object
        * @param axis Axis of the histogram
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_correlationStyle(TH1D* hist, const std::string& axis, double range, int rebin);

        /**
        * @brief Set charge distribution style
        * @param hist Pointer to TH1D object
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_chargeStyle(TH1D* hist, double range, int rebin);

        /**
        * @brief Set residual distribution style
        * @param hist Pointer to TH1D object
        * @param axis Axis of the histogram
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_residualStyle(TH1D* hist, const std::string& axis, double range, int rebin);

        /**
        * @brief Set track chi2 distribution style
        * @param hist Pointer to TH1D object
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_trackChi2Style(TH1D* hist, double range, int rebin);

        /**
        * @brief Set track chi2/ndof distribution style
        * @param hist Pointer to TH1D object
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_trackChi2ndofStyle(TH1D* hist, double range, int rebin);

        /**
        * @brief Set track angle distribution style
        * @param hist Pointer to TH1D object
        * @param axis Axis of the histogram
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_trackAngleStyle(TH1D* hist, const std::string& axis, double range, int rebin);

        /**
        * @brief Set cluster multiplicity per event distribution style
        * @param hist Pointer to TH1D object
        * @param range Range of the histogram
        * @param rebin Rebinning factor
        **/
        void set_multiplicityStyle(TH1D* hist, double range, int rebin);

        /**
        * @brief Set hitmap projection style
        * @param hist Pointer to TH1D object
        * @param axis Axis of the histogram
        * @param rebin Rebinning factor
        **/
        void set_projectionStyle(TH1D* hist, const std::string& axis, int rebin);

        /**
        * @brief Set hitmap projection style
        * @param hist Pointer to TH1D object
        * @param axis Axis of the histogram
        * @param rebin Rebinning factor
        **/
        void write_allRefProjection(TFile* input_file,
                                    TFile* output_file,
                                    std::vector<std::string> refnames,
                                    std::string axis);

        std::vector<std::string> refnames, dutnames;
        std::string input_filename, output_filename;
        std::vector<std::string> rootFile_dir_names;
        std::vector<TDirectory*> rootFile_dirs;
        std::string input_fileDir, base_filename;
        TFile *input_file, *output_file;
        TDirectory *hitmap, *correlation, *analysisRef, *track, *analysisCE65;
};