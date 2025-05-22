#ifndef TRACK_RESOLUTION_HPP
#define TRACK_RESOLUTION_HPP

#include "plot_histogram/plot_histogram.hpp"

struct detector_information {
    std::string detector_name;
    std::string detector_type;
    std::string detector_role;
    std::vector<double> xyz_position;
    std::vector<double> xyz_orientation;
    std::vector<int> xy_pixels;
    std::vector<double> position_resolution;
    double material_budget;
    double time_resolution;
};

class track_resolution : public plot_histogram {
    public:
        /**
        * @brief Constructor for this track_resolution object
        **/
        track_resolution();

        // function
        /**
        * @brief Calculate mean value of the data
        * @param data Vector of data
        **/
        double calculate_mean(std::vector<double> data);

        double calculate_extrapolationResolution(const std::vector<double>& ref_position, const std::vector<double>& ref_resolution, double extr_position);

        /**
        * @brief calculate track resolution
        * @param ref_resolution Vector of reference resolution
        * @param ref_position Vector of reference position
        * @param dut_position Position of DUT
        **/
        double calculate_trackResolution(std::vector<double> ref_resolution, std::vector<double> ref_position, double dut_position);

        /**
        * @brief Open .conf file
        * @param filename Name of the configuration file
        **/
        std::ifstream open_confFile(std::string filename);

        /**
        **/
        std::vector<detector_information> import_detectorInformation(std::ifstream& file);

        void run();

    private:
        // parameter
        
};

#endif // TRACK_RESOLUTION_HPP