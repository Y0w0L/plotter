#include "track_resolution.hpp"

track_resolution::track_resolution() {
    std::cout << "track_resolution object is created" << std::endl;
}

double track_resolution::calculate_mean(std::vector<double> data) {
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

double track_resolution::calculate_extrapolationResolution(const std::vector<double>& ref_position, const std::vector<double>& ref_resolution, double extr_position) {
    std::cout << "Calculate extrapolation resolution" << std::endl;
    int ref_size = ref_position.size();
    double sum_weight = 0;
    double sum_weight_position = 0;
    double sum_weight_position2 = 0;
    
    for(int i = 0; i < ref_size; i++) {
        double weight = 1 / pow(ref_resolution[i], 2);
        sum_weight += weight;
        sum_weight_position += weight * ref_position[i];
        sum_weight_position2 += weight * pow(ref_position[i], 2);
        std::cout << "Weight: " << weight << std::endl;
    }
    
    double delta = sum_weight * sum_weight_position2 - pow(sum_weight_position, 2);
    double sigma_a2 = sum_weight_position2 / delta;
    double sigma_b2 = sum_weight / delta;
    double cov_ab = -sum_weight_position / delta;
    std::cout << "Delta: " << delta << std::endl;
    std::cout << "Sigma_a2: " << sigma_a2 << std::endl;
    std::cout << "Sigma_b2: " << sigma_b2 << std::endl;
    std::cout << "Cov_ab: " << cov_ab << std::endl;

    double extr_resolution = sqrt(sigma_b2 + pow(extr_position, 2) * sigma_a2 + 2 * extr_position * cov_ab);
    std::cout << "2: " << pow(extr_position, 2) * sigma_a2 << std::endl;
    std::cout << "3: " << 2 * extr_position * cov_ab << std::endl;
    return extr_resolution;
}

double track_resolution::calculate_trackResolution(std::vector<double> ref_resolution, std::vector<double> ref_position, double dut_position) {
    double mean_position;
    double mean_resolution;
    double exp_mean_position;
    double variance_position = 0;
    double track_resolution;

    mean_position = calculate_mean(ref_position);
    mean_resolution = calculate_mean(ref_resolution);
    exp_mean_position = dut_position - mean_position;
    for(std::vector<double>::iterator it = ref_position.begin(); it != ref_position.end(); it++) {
        variance_position += pow(*it - mean_position, 2);
    }
    track_resolution = mean_resolution * sqrt(1 + pow(exp_mean_position, 2) / variance_position);
    std::cout << "Mean position: " << mean_position << std::endl;
    std::cout << "Expected mean position: " << exp_mean_position << std::endl;
    std::cout << "Variance position: " << variance_position << std::endl;
    return track_resolution;
}

// Open file
std::ifstream track_resolution::open_confFile(std::string filename) {
    std::ifstream file(filename);
    std::cout << "Open file: " << filename << std::endl;
    if(!file.is_open()) {
        print_message(Form("Error: Cannot open config file in %s", filename.c_str()), RED);
    }
    return std::move(file);
}

std::vector<detector_information> track_resolution::import_detectorInformation(std::ifstream& file) {
    std::cout << "Import detector information" << std::endl;
    std::vector<detector_information> detector_info;
    std::string line;
    std::string current_detectorName;
    bool isDUT = false;
    bool isREF = false;

    while(std::getline(file, line)) {
        if(line.empty() || line[0] == '#') {
            current_detectorName = line.substr(1, line.size() - 2);
            isDUT = (current_detectorName.find("CE65") != std::string::npos);
            isREF = (current_detectorName.find("ALPIDE") != std::string::npos);
            continue;
        }

        std::istringstream iss(line);
        std::string key, equal, value;
        if(iss >> key >> equal >> value) {}
    }

    return detector_info;
}

void track_resolution::run() {
    // std::ifstream file = track.open_confFile("detector.conf");
    // track.import_detectorInformation(file);

    //example
    std::vector<double> ref_resolution = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006}; // um
    std::vector<double> ref_position = {0, 25.4, 50.8, 127, 152.4, 177.8}; // mm
    double dut_position = 76.2;

    double track_resolution = track_resolution::calculate_trackResolution(ref_resolution, ref_position, dut_position);
    std::cout << "Track resolution: " << track_resolution << std::endl;

    //std::vector<double> ref_resolutions = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006}; // um
    std::vector<double> ref_resolutions(6, 0.006);
    double extr_resolution = track_resolution::calculate_extrapolationResolution(ref_position, ref_resolutions, dut_position);
    std::cout << "Extrapolation resolution: " << extr_resolution * 1000 << std::endl;
}