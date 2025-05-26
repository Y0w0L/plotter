#include "track_resolution.h"

track_resolution::track_resolution() {
    std::cout << "track_resolution object is created" << std::endl;
}

double track_resolution::calculate_mean(const std::vector<double>& data) {
    if(data.empty()) {
        LOG_ERROR.source("track_resolution::calculate_mean") << "Input data vector is empty.";
        throw std::invalid_argument("Input data vector cannot be empty for calculata_mean.");
    }
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

// calculation of track resolution
// double calculate_straightlineResolution(std::vector<double> ref_resolution, std::vector<double> ref_position) {
//     double mean_position;
//     double mean_resolution;
//     double ref_mean_position;

//     mean_position = calculate_mean(ref_position);
//     mean_resolution = calculate_mean(ref_resolution);
//     while(std::vector<double> ::iterator it = ref_position.begin(); it != ref_position.end(); it++) {
//         ref_mean_position += *it - mean_position;
//     }
// }

double track_resolution::calculate_extrapolationResolution(const std::vector<double>& ref_position,
                                                           const std::vector<double>& ref_resolution,
                                                           double extr_position) {
    LOG_DEBUG.source("track_resolution::calculate_extrapolationResolution") << "Calculate extrapolation resolution.";
    if(ref_position.size() != ref_resolution.size() || ref_resolution.empty()) {
        throw std::invalid_argument("Input vectors mismatch or empty for calculate_extrapolationResolution.");
    }
    
    int ref_size = ref_position.size();
    double sum_weight = 0;
    double sum_weight_position = 0;
    double sum_weight_position2 = 0;
    double weight;
    
    for(int i = 0; i < ref_size; i++) {
        if(ref_resolution[i] == 0.0) {
            throw std::runtime_error("Division by zero: ref_resolution cannot be zero.");
        }
        weight = 0.0;
        weight = 1 / pow(ref_resolution[i], 2);
        sum_weight += weight;
        sum_weight_position += weight * ref_position[i];
        sum_weight_position2 += weight * pow(ref_position[i], 2);
        LOG_DEBUG.source("track_resolution::calculate_extrapolationResolution") << "Weight: " << weight;
    }
    
    double delta = sum_weight * sum_weight_position2 - pow(sum_weight_position, 2);
    if(delta == 0) {
        throw std::runtime_error("Divison by zero: Delta is zero in calculate_extrapolationResolution.");
    }
    
    double sigma_a2_variance_intercept = sum_weight_position2 / delta; // Renamed for clarity
    double sigma_b2_variance_slope = sum_weight / delta;            // Renamed for clarity
    double cov_ab_covariance_slope_intercept = -sum_weight_position / delta; // Renamed for clarity
    LOG_DEBUG.source("track_resolution::calculate_extrapolationResolution") << "Delta: " << delta;
    LOG_DEBUG.source("track_resolution::calculate_extrapolationResolution") << "Sigma_a2: " << sigma_a2_variance_intercept;
    LOG_DEBUG.source("track_resolution::calculate_extrapolationResolution") << "Sigma_b2: " << sigma_b2_variance_slope;
    LOG_DEBUG.source("track_resolution::calculate_extrapolationResolution") << "Cov_ab: " << cov_ab_covariance_slope_intercept; 

    double term1_slope_variance_contribution = sigma_b2_variance_slope;
    double term2_intercept_variance_contribution = (extr_position * extr_position) * sigma_a2_variance_intercept; // Replaced pow(x,2)
    double term3_covariance_contribution = 2.0 * extr_position * cov_ab_covariance_slope_intercept;

    double extr_resolution_squared = term1_slope_variance_contribution + term2_intercept_variance_contribution + term3_covariance_contribution;
    if(extr_resolution_squared <  0) {
        LOG_WARNING.source("track_resolution::calculate_extrapolationResolution") << "Negative value under square root for extrapolation resolution. Result will be NaN.";
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(extr_resolution_squared);
}

double track_resolution::calculate_trackResolution(const std::vector<double>& ref_resolution,
                                                   const std::vector<double>& ref_position,
                                                   double dut_position) {
    if(ref_position.empty() || ref_resolution.empty() || ref_position.size() != ref_resolution.size()) {
        throw std::invalid_argument("Input vectors mismatch or empty for calculate_trackResolution.");
    }
                                                    
    double mean_position = calculate_mean(ref_position);
    double mean_resolution = calculate_mean(ref_resolution);
    double exp_mean_position_diff = dut_position - mean_position;
    double diff;
    double variance_position = 0;

    for(double pos : ref_position) {
        diff = 0.0;
        diff = pos - mean_position;
        variance_position += diff * diff;
    }
    if(variance_position == 0.0) {
        throw std::runtime_error("Division by zero: variance_position is zero in calculate_trackResolution.");
    }

    double term_inside_sqrt = 1.0 + (exp_mean_position_diff * exp_mean_position_diff) / variance_position;
    if(term_inside_sqrt < 0) {
        LOG_WARNING.source("track_resolution::calculate_trackResolution") << "Negative value under square root for track resolution calculation. Result will be NaN.";
        return std::numeric_limits<double>::quiet_NaN();
    }

    double track_res = mean_resolution * std::sqrt(term_inside_sqrt);
    LOG_DEBUG.source("track_resolution::calculate_trackResolution") << "Mean position: " << mean_position;
    //LOG_DEBUG.source("track_resolution::calculate_trackResolution") << "Expected mean position: " >> exp_mean_position_diff;
    LOG_DEBUG.source("track_resolution::calculate_trackResolution") << "Variance position: " << variance_position;
    LOG_DEBUG.source("track_resolution::calculate_trackResolution") << "Track resolution: " << track_res;

    return track_res;
}

// Open file
std::ifstream track_resolution::open_confFile(const std::string& filename) {
    LOG_DEBUG.source("track_resolution::open_confFile") << "Opening gile: " << filename;
    std::ifstream file(filename);
    if(!file.is_open()) {
        LOG_FATAL.source("track_resolution::open_confFile") << "Cannot open config file: " << filename;
    }
    return file;
}

// std::vector<detector_information> track_resolution::import_detectorInformation(std::ifstream& file) {
//     std::cout << "Import detector information" << std::endl;
//     std::vector<detector_information> detector_info;
//     std::string line;
//     std::string current_detectorName;
//     bool isDUT = false;
//     bool isREF = false;

//     while(std::getline(file, line)) {
//         if(line.empty() || line[0] == '#') {
//             current_detectorName = line.substr(1, line.size() - 2);
//             isDUT = (current_detectorName.find("CE65") != std::string::npos);
//             isREF = (current_detectorName.find("ALPIDE") != std::string::npos);
//             continue;
//         }

//         std::istringstream iss(line);
//         std::string key, equal, value;
//         if(iss >> key >> equal >> value) {}
//     }

//     return detector_info;
// }

void track_resolution::run() {
    // std::ifstream file = track.open_confFile("detector.conf");
    // track.import_detectorInformation(file);

    //example
    std::vector<double> ref_resolution = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006}; // um
    std::vector<double> ref_position = {0, 25.4, 50.8, 127, 152.4, 177.8}; // mm
    double dut_position = 76.2;

    double track_resolution = track_resolution::calculate_trackResolution(ref_resolution, ref_position, dut_position);
    LOG_INFO.source("track_resolution::run") << "Track resolution: " << track_resolution << "um";
    //std::vector<double> ref_resolutions = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006}; // um
    std::vector<double> ref_resolutions(6, 0.006);
    double extr_resolution = track_resolution::calculate_extrapolationResolution(ref_position, ref_resolutions, dut_position);
    //std::cout << "Extrapolation resolution: " << extr_resolution * 1000 << std::endl;
    LOG_INFO.source("track_resolution::run") << "Extrapolation resolution: " << extr_resolution * 1000;
}