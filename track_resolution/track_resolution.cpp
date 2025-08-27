#include "track_resolution.h"

const double RAD_LENGTH_SILICON = 93.7; //mm
const double RAD_LENGTH_AIR = 304200.0;; // mm STP
const double RAD_LENGTH_ALUMINUM = 89.0; //mm
const double RAD_LENGTH_MYLAR = 286.0; // mm PET
const double RAD_LENGTH_SCINTI = 438.0; // mm

const Particle ELECTRON = {"electron", 0.511, -1};
const Particle PION_M = {"pi-", 139.6, -1};
const Particle PROTON = {"proton", 938.3, 1};

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

double track_resolution::calculate_beta(
    double momentum_MeV_c,
    double mass_MeV_c2
) {
    if(momentum_MeV_c <= 0) return 0.0;
    double total_energy_MeV = std::sqrt(momentum_MeV_c*momentum_MeV_c + mass_MeV_c2*mass_MeV_c2);
    if(total_energy_MeV == 0) return 0.0;

    return momentum_MeV_c / total_energy_MeV;
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

double track_resolution::calculate_extrapolationResolution(
    const std::vector<double>& ref_position,
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

    // double term1_slope_variance_contribution = sigma_b2_variance_slope;
    // double term2_intercept_variance_contribution = (extr_position * extr_position) * sigma_a2_variance_intercept; // Replaced pow(x,2)
    // double term3_covariance_contribution = 2.0 * extr_position * cov_ab_covariance_slope_intercept;

    //double extr_resolution_squared = term1_slope_variance_contribution + term2_intercept_variance_contribution + term3_covariance_contribution;
    double extr_resolution_squared = (extr_position * extr_position) * sigma_b2_variance_slope + sigma_a2_variance_intercept + 2.0 * extr_position * cov_ab_covariance_slope_intercept;
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

    const size_t ref_size = ref_position.size();

    for(double pos : ref_position) {
        diff = 0.0;
        diff = pos - mean_position;
        variance_position += diff * diff;
    }
    if(variance_position == 0.0) {
        throw std::runtime_error("Division by zero: variance_position is zero in calculate_trackResolution.");
    }

    double term_inside_sqrt = (1.0 / ref_size) + (exp_mean_position_diff * exp_mean_position_diff) / variance_position;
    //double term_inside_sqrt = (1.0 / ref_size) + variance_position / (exp_mean_position_diff * exp_mean_position_diff);
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

void track_resolution::calculate_tResoltuion(const std::vector<double>& ref_resolution, const std::vector<double>& ref_position, double dut_position) {
    LOG_DEBUG.source("track_resolution::calculate_tResolution") << "Start calculation of track resolution.";
    double sum_0 = 0.0;
    double sum_1 = 0.0;
    double sum_2 = 0.0;
    double D = 0.0;
    double track_resolution = 0.0;
    
    double temp_value;
    for(int i=0; i< ref_position.size(); i++) {
        temp_value = 1 / (ref_resolution[i]*ref_resolution[i]);
        sum_0 += temp_value;
        sum_1 += ref_position[i] * temp_value;
        sum_2 -= ref_position[i]*ref_position[i] * temp_value;
    }
    D = sum_0*sum_2 - sum_1*sum_1;
    track_resolution = std::sqrt(1/D * (sum_2 + dut_position*dut_position*sum_0 - 2*dut_position*sum_1));
    std::cout << "Track resolution-> " << track_resolution << std::endl;
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

// calculate theta_0 by multiple scattering
double track_resolution::calculate_theta0(double momentum, double beta, int charge_z, double thickness, double radiation_length) {
        if (momentum <= 0 || beta <= 0 || thickness <= 0 || radiation_length <= 0) {
        return 0.0;
    }
    double x_over_X0 = thickness / radiation_length;
    if (x_over_X0 <= 0) return 0.0;
    
    double log_term = 1.0 + 0.038 * std::log(x_over_X0);
    return (13.6 / (beta * momentum)) * charge_z * std::sqrt(x_over_X0) * log_term;
}

std::vector<ScatteringEffect> track_resolution::calculate_multilayer_scattering(
    double momentum, double beta, int charge_z,
    const std::vector<MaterialLayer>& layers) {
    
    std::vector<ScatteringEffect> results_per_layer;
    ScatteringConvariance sigma = {0.0, 0.0, 0.0};
    double theta_i, L, theta_i2;
    ScatteringEffect current_effect;

    for(const auto& layer : layers) {
        L = layer.thickness_mm;

        sigma.yy = sigma.yy + 2*L*sigma.yth + L*L*sigma.thth;
        sigma.yth = sigma.yth + L*sigma.thth;
        // sigma.thth is constant
        theta_i = calculate_theta0(momentum, beta, charge_z, L, layer.rad_length_mm);

        theta_i2 = theta_i*theta_i;
        sigma.yy += (L*L/3.0) * (theta_i2);
        sigma.yth += (L/2.0) * (theta_i2);
        sigma.thth += (theta_i2);

        current_effect.angle_rms_rad = std::sqrt(sigma.thth);
        current_effect.displacement_rms_mm = std::sqrt(sigma.yy);
        results_per_layer.push_back(current_effect);
    }

    return results_per_layer;
    }

void track_resolution::get_multipleScatteringEffect() {
    LOG_STATUS.source("track_resolution::get_multipleScatteringEffect") << "Start get_multipleScatteringEffect.";
    std::vector<double> distances = {0.01, 0.05, 0.1, 1, 10, 25, 50, 100, 200, 500, 1000}; // mm
    std::map<std::string, double> radiationLengths = {
        {"Air", RAD_LENGTH_AIR},
        {"Scinti", RAD_LENGTH_SCINTI},
        {"Mylar", RAD_LENGTH_MYLAR},
        {"Silicon", RAD_LENGTH_SILICON},
        {"Aluminum", RAD_LENGTH_ALUMINUM}
    };

    std::vector<std::string> materialOrder = {"Air", "Scinti", "Mylar", "Silicon", "Aluminum"};

    const double particle_momentum =120000; // Mev/c
    const double beta = 1;
    const int z = 1;

    std::map<std::string, TGraph*> angle_graphs;
    std::map<std::string, TGraph*> distance_graphs;
    std::map<double, std::map<std::string, double>> angle_results;
    std::map<double, std::map<std::string, double>> distance_results;

    for(const auto& material : radiationLengths) {
        const std::string& name = material.first;
        angle_graphs[name] = new TGraph();
        distance_graphs[name] = new TGraph();
    }

    for(const auto& material : radiationLengths) {
        const std::string& name = material.first;
        const double radiation_length = material.second;

        // angle_graphs[name] = new TGraph;
        // distance_graphs[name] = new TGraph;
        int pointIndex = 0;
        for(const double dist : distances) {
            double angle = track_resolution::calculate_theta0(particle_momentum, beta, z, dist, radiation_length);
            double distance = std::sqrt((dist*dist/3.0) * (angle*angle));
            angle_graphs[name]->SetPoint(pointIndex++, dist, angle);
            distance_graphs[name]->SetPoint(pointIndex++, dist, distance);
            angle_results[dist][name] = angle;
            distance_results[dist][name] = distance;
        }
    }

    TCanvas* canvas = new TCanvas("c", "c", 800, 600);
    // canvas->SetLogx();
    canvas->SetLogy();
    canvas->SetGrid();
    canvas->SetTopMargin(0.062);
    canvas->SetBottomMargin(0.14);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.03);

    TMultiGraph* angle_mg = new TMultiGraph();
    angle_mg->SetTitle(";distance [mm];scattering angle [rad]");
    std::map<std::string, int> colors = {
        {"Air", kGray + 1}, {"Scintillator", kGreen + 1}, {"Mylar", kBlue},
        {"Silicon", kOrange + 7}, {"Aluminum", kRed}
    };
    std::map<std::string, int> markers = {
        {"Air", 20}, {"Scintillator", 21}, {"Mylar", 22},
        {"Silicon", 23}, {"Aluminum", 34}
    };
    TLegend* legend = new TLegend(0.15, 0.65, 0.45, 0.88);
    legend->SetHeader("Materials");
    legend->SetTextSize(0.025);

    for (const auto& name : materialOrder) {
        TGraph* g = angle_graphs[name];
        g->SetLineColor(colors[name]);
        g->SetLineWidth(2);
        g->SetMarkerColor(colors[name]);
        g->SetMarkerStyle(markers[name]);
        g->SetMarkerSize(1.2);
        
        angle_mg->Add(g, "P");
        legend->AddEntry(g, name.c_str(), "lp");
    }

    // angle_mg->GetXaxis()->SetRangeUser(0.01, 1000);
    // angle_mg->GetYaxis()->SetRangeUser(0.001, )
    angle_mg->Draw("A");
    canvas->SaveAs("multiple_scattering_angle.pdf");
}

void track_resolution::run() {
    // std::ifstream file = track.open_confFile("detector.conf");
    // track.import_detectorInformation(file);

    //example
    std::vector<double> ref_resolution(6, 0.007);
    std::vector<double> ref_position;
    // ref_resolution = {0.007, 0.007, 0.007, 0.007, 0.007, 0.007}; // mm
    //ref_position = {0, 25.4, 50.8, 127, 152.4, 177.8}; // mm
    ref_position = {0, 25, 50, 150, 175, 200}; // mm
    //ref_position = {0, 25.4, 50.8, 127, 152.4}; // mm
    //ref_position = {25.4, 50.8, 127, 152.4};
    //double dut_position = 76.2;
    double dut_position = 100;

    double track_resolution = track_resolution::calculate_trackResolution(ref_resolution, ref_position, dut_position);
    LOG_INFO.source("track_resolution::run") << "Track resolution: " << track_resolution * 1000 << " um";
    //std::vector<double> ref_resolutions = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006}; // um
    // std::vector<double> ref_resolutions(6, 0.006);
    double extr_resolution = track_resolution::calculate_extrapolationResolution(ref_position, ref_resolution, dut_position);
    //std::cout << "Extrapolation resolution: " << extr_resolution * 1000 << std::endl;
    LOG_INFO.source("track_resolution::run") << "Extrapolation resolution: " << extr_resolution * 1000 << " um";

    LOG_STATUS.source("track_resolution::run") << "Start considering multiple scattering effect.";

    //double p = 3000.0; //momentum 3 GeV/c
    double p = 120000.0;
    double beta = 1.0; // almost speed of light
    int z = 1; // charge

    const double mass_electron = 0.511;
    const double mass_pion = 140;
    const double mass_proton = 938.3;

    auto calculate_beta = [&](double mass) {
        return p / std::sqrt(p*p + mass*mass);
    };

    double beta_electron  = calculate_beta(mass_electron);
    double beta_pion      = calculate_beta(mass_pion);
    double beta_proton    = calculate_beta(mass_proton);
    std::cout << beta_pion << std::endl;

    std::vector<MaterialLayer> geometry = {
        {"ALPIDE_0",   0.05, RAD_LENGTH_SILICON}, // Silicon: 50um
        {"Air_0",       25.4, RAD_LENGTH_AIR}, // Air: 2.54cm
        {"ALPIDE_1",   0.05, RAD_LENGTH_SILICON}, // Silicon: 50um
        {"Air_1",       25.4, RAD_LENGTH_AIR},
        {"ALPIDE_2",   0.05, RAD_LENGTH_SILICON},
        {"Air_2",       25.4, RAD_LENGTH_AIR},
        {"CE65_3",   0.05, RAD_LENGTH_SILICON}, // Silicon: 50um
        {"Air_3",       25.4, RAD_LENGTH_AIR}, // Air: 2.54cm
        {"DPTS",   0.05, RAD_LENGTH_SILICON}, // Silicon: 50um
        {"Air_4",       25.4, RAD_LENGTH_AIR},
        {"ALPIDE_4",   0.05, RAD_LENGTH_SILICON},
        {"Air_5",       25.4, RAD_LENGTH_AIR}, // Air: 2.54cm
        {"ALPIDE_5",   0.05, RAD_LENGTH_SILICON}, // Silicon: 50um
        {"Air_6",       25.4, RAD_LENGTH_AIR},
        {"ALPIDE_6",   0.05, RAD_LENGTH_SILICON},
    };

    //std::vector<ScatteringEffect> effects = this->calculate_multilayer_scattering(p, beta, z, geometry);
    std::vector<ScatteringEffect> effects_electron = this->calculate_multilayer_scattering(p, beta_electron, z, geometry);
    std::vector<ScatteringEffect> effects_pion = this->calculate_multilayer_scattering(p, beta_pion, z, geometry);
    std::vector<ScatteringEffect> effects_proton = this->calculate_multilayer_scattering(p, beta_proton, z, geometry);


    LOG_STATUS.source("track_resolution::run") << p / 1000 << " GeV/c, z = " << z;
    LOG_STATUS.source("track_resolution::run") << "=====================================================================";
    for(size_t i=0; i<geometry.size(); i++) {
        LOG_STATUS.source("track_resolution::run") << " - " << geometry[i].name << " Thickness: " << geometry[i].thickness_mm << " mm";
        LOG_STATUS.source("track_resolution::run") << ">>> Angle RMS       : " << effects_pion[i].angle_rms_rad * 1000 << " mrad";
        LOG_STATUS.source("track_resolution::run") << ">>> Displacement RMS: " << effects_pion[i].displacement_rms_mm * 1000.0 << " um";
        LOG_STATUS.source("track_resolution::run") << "-----------------------------------------------------------------";
    }
    LOG_STATUS.source("track_resolution::run") << "Track resolution: " << track_resolution * 1000.0 << " um";
    LOG_STATUS.source("track_resolution::run") << "Multiple scattering effect: " << effects_pion[5].displacement_rms_mm * 1000.0 << " um";
    double TR2 = (track_resolution * 1000.0) * (track_resolution * 1000.0);
    double MSE2 = (effects_pion[5].displacement_rms_mm * 1000.0)*(effects_pion[5].displacement_rms_mm * 1000.0);
    LOG_STATUS.source("track_resolution::run") << "TR^2: " << TR2;
    LOG_STATUS.source("track_resolution::run") << "MSE^2: " << MSE2;
    LOG_STATUS.source("track_resolution::run") << "SUM(TR^2 + MSE^2): " << TR2 + MSE2;

    // get_multipleScatteringEffect();
    track_resolution::calculate_tResoltuion(ref_resolution, ref_position, dut_position);
}   