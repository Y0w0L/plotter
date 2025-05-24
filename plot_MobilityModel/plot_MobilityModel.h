#ifndef PLOT_MOBILITYMODEL_H
#define PLOT_MOBILITYMODEL_H

#include "plot_histogram/plot_histogram.h"

class plot_MobilityModel : public plot_histogram {
    public:
        /**
        * @brief Constructor for this plot_MobilityModel object
        **/
        plot_MobilityModel();

        // function
        /**
        * @brief Set cluster size style
        * @param hist Pointer to TH1D object
        * @param color Marker color
        * @param marker_style Marker style
        **/
        static void clSizeStyle(TH1D* hist, const int color, const int marker_style);

        /**
        * @brief Get histogram from TFile object
        * @param file Pointer to TFile object
        * @param hist_name Name of the histogram
        **/
        static TH1D* get_hist(TFile* file, const std::string& hist_name);

        /**
        * @brief Get TFile objects from the list of filenames
        * @param filenames List of filenames
        **/
        static std::vector<TFile*> get_TFiles(const std::vector<std::string>& filenames);

        /**
        * @brief Get TH1D objects from the list of TFile objects
        * @param files List of TFile objects
        * @param hist_name Name of the histogram
        **/
        static std::vector<TH1D*> get_hists(std::vector<TFile*>& files, const std::string& hist_name);

        /**
        * @brief Get fits from the list of TH1D objects
        * @param hists List of TH1D objects
        **/
        static std::vector<TF1*> get_fits(std::vector<TH1D*> hists);

        /**
        * @brief Get Landau fits from the list of TH1D objects
        * @param hists List of TH1D objects
        **/
        static std::vector<TF1*> get_landauFits(std::vector<TH1D*> hists);

        /**
        * @brief Write cluster size plots
        * @param hists List of TH1D objects
        * @param names List of names
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        **/
        static void write_clSizePlots(std::vector<TH1D*> hists, 
                                      const std::vector<std::string>& names,
                                      const int& x_max,
                                      const int& y_max);

        /**
        * @brief Write charge plots
        * @param hists List of TH1D objects
        * @param names List of names
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        **/
        static std::vector<TH1D*> write_chargePlots(std::vector<TH1D*> hists,
                                                    std::vector<std::string> names,
                                                    int x_max,
                                                    int y_max);

        /**
        * @brief Write residual plots
        * @param hists List of TH1D objects
        * @param fits List of TF1 objects
        * @param names List of names
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        * @param axis Axis name
        **/
        static void write_residualPlots(std::vector<TH1D*> hists,
                                        std::vector<TF1*> fits,
                                        std::vector<std::string> names,
                                        int x_max,
                                        int y_max,
                                        std::string axis);

        /**
        * @brief Write in-pixel residual plots
        * @param hists List of TH1D objects
        * @param names List of names
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        * @param axis Axis name
        **/
        static void write_inPixelResidualPlots(std::vector<TH1D*> hists,
                                               std::vector<std::string> names,
                                               int x_max,
                                               int y_max,
                                               std::string axis);

        /**
        * @brief Save cluster size plots
        * @param output_file Pointer to the output ROOT file
        * @param filenames List of filenames
        * @param names List of names
        * @param save_name Name of the saved file
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        **/
        static void save_clSizePlots(TFile* output_file,
                                     std::vector<std::string> filenames,
                                     std::vector<std::string> names,
                                     std::string save_name,
                                     int x_max,
                                     int y_max);
        
        /**
        * @brief Save charge plots
        * @param output_file Pointer to the output ROOT file
        * @param filenames List of filenames
        * @param names List of names
        * @param save_name Name of the saved file
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        **/
        static void save_chargePlots(TFile* output_file,
                                     std::vector<std::string> filenames,
                                     std::vector<std::string> names,
                                     std::string save_name,
                                     int x_max,
                                     int y_max);

        static void save_seedChargePlots(TFile* output_file,
                                         std::vector<std::string> filenames,
                                         std::vector<std::string> names,
                                         std::string save_name,
                                         int x_max,
                                         int y_max);
        
        /**
        * @brief Save residual plots
        * @param output_file Pointer to the output ROOT file
        * @param filenames List of filenames
        * @param names List of names
        * @param save_name Name of the saved file
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        **/
        static void save_residualPlots(TFile* output_file,
                                       std::vector<std::string> filenames,
                                       std::vector<std::string> names,
                                       std::string save_name,
                                       int x_max,
                                       int y_max);
        
        /**
        * @brief Save in-pixel residual plots
        * @param output_file Pointer to the output ROOT file
        * @param filenames List of filenames
        * @param names List of names
        * @param save_name Name of the saved file
        * @param x_max Maximum value of x-axis
        * @param y_max Maximum value of y-axis
        **/
        static void save_inPixelResidualPlots(TFile* output_file,
                                              std::vector<std::string> filenames,
                                              std::vector<std::string> names,
                                              std::string save_name,
                                              int x_max,
                                              int y_max);

        static void save_suminPixelResidualPlots(TFile* output_file,
                                            std::vector<std::string> filenames,
                                            std::vector<std::string> names,
                                            std::string save_name,
                                            int x_max,
                                            int y_max);

        void run();

    private:
        // parameter
        std::string input_fileDir, output_fileDir, output_file;

        std::vector<std::string> pixel_pitchs;
        std::vector<std::string> chip_types;
        std::vector<std::string> chip_variations;
        std::vector<std::string> threshold_values;
        std::vector<std::string> thresholds;
        std::vector<std::string> model_names;
        std::vector<std::string> models;
        
        TH1D* h_clSize;
        int clSize_x_max, clSize_x_max_highThd, charge_x_max;
        int clSize_y_max, clSize_y_max_highThd, charge_y_max;
        int x_residual_max, y_residual_max;
        int y_inPixelResidual_max;
};

#endif