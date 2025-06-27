#ifndef TRACK_RESOLUTION_HPP
#define TRACK_RESOLUTION_HPP

#include "plot_histogram/plot_histogram.h"
#include "Messenger/Messenger.h"

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

struct MaterialLayer {
    std::string name;
    double thickness_mm;
    double rad_length_mm;
};

struct ScatteringConvariance {
    double yy;
    double yth;
    double thth;
};

struct ScatteringEffect {
    double angle_rms_rad;
    double displacement_rms_mm;
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
        double calculate_mean(const std::vector<double>& data);

        double calculate_extrapolationResolution(const std::vector<double>& ref_position, const std::vector<double>& ref_resolution, double extr_position);

        /**
        * @brief calculate track resolution
        * @param ref_resolution Vector of reference resolution
        * @param ref_position Vector of reference position
        * @param dut_position Position of DUT
        **/
        double calculate_trackResolution(const std::vector<double>& ref_resolution, const std::vector<double>& ref_position, double dut_position);

        void calculate_tResoltuion(const std::vector<double>& ref_resolution, const std::vector<double>& ref_position, double dut_position);

        /**
        * @brief Open .conf file
        * @param filename Name of the configuration file
        **/
        std::ifstream open_confFile(const std::string& filename);

        double calculate_theta0(double momentum, double beta, int charge_z, double thickness, double radiation_length);
        std::vector<ScatteringEffect> calculate_multilayer_scattering(double momentum, double beta, int charge_z, const std::vector<MaterialLayer>& layers);
        void get_multipleScatteringEffect();
        /**
        **/
        // std::vector<detector_information> import_detectorInformation(std::ifstream& file);

        void run();

    private:
        // parameter
        
};

#endif // TRACK_RESOLUTION_HPP