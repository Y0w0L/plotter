// #include "plot_AllpixData.h"

// using namespace allpix;

// plot_AllpixData::plot_AllpixData() {
//     LOG_STATUS.source("plot_AllpixData::plot_AllpixData") << "plot_AllpixData object is created.";
//     gSystem->Load("/home/towa/package/allpix/install/lib/libAllpixObjects.so");
// }

// plot_AllpixData::~plot_AllpixData() {

// }

// void plot_AllpixData::run_analysis() {
//     // if(config.empty()) {
//     //     LOG_WARNING.source("plot_AllpixData::run_analysis") << "No chip configuration provided.";
//     //     return;
//     // }

//     // h_residual_x = new TH1D("residual_x", "residual_x", 200, -50, 50);
//     // h_residual_y = new TH1D("residual_y", "residual_y", 200, -50, 50);
//     // h_residual_x_inPixel = new TH1D("residual_x_inPixel", "residual_x_inPixel",
//     //                                 20, -config.pixel_pitch/2, config.pixel_pitch/2,
//     //                                 20, -config.pixel_pitch/2, config.pixel_pitch/2);

//     // std::string input_file_name = "/home/towa/alice3/data/ce65_sim_202505/n10v/ce65_p15_gap_Thd0e_e3GeV_masetti.root";

//     // const double seed_threshold = 400;
//     // const double neighbor_threshold = 10;

//     // const std::string mc_branch_name = "MCParticle";
//     // const std::string pixel_branch_name = "PixelHit";

//     // TFile* input_file = TFile::Open(input_file_name.c_str());
//     // if(!input_file || input_file->IsZombie()) {
//     //     LOG_ERROR.source("plot_AllpixData::run_analysis") << "Cannot open input file " << input_file_name;
//     //     return;
//     // }

//     // TTree* tree = dynamic_cast<TTree*>(input_file->Get("tree"));
//     // if(!tree) {
//     //     LOG_ERROR.source("plot_AllpixData::run_analysis") << "Cannot find TTree 'tree' in the file.";
//     //     return; 
//     // }

//     // std::vector<allpix::MCParticle>* mc_particles = nullptr;
//     // std::vector<allpix::PixelHit>* pixel_hits = nullptr;

//     // tree->SetBranchAddress(mc_branch_name.c_str(), &mc_particles);
//     // tree->SetBranchAddress(pixel_branch_name.c_str(), &pixel_hits);


    
// }

// std::vector<Cluster> plot_AllpixData::doClustering(const std::vector<PixelHit*>& pixel_hits) const {
//     std::vector<Cluster> clusters;
//     std::map<const PixelHit*, bool> usedPixel;

//     // std::vector<const PixelHit*> pixel_hits;
//     // for(const auto& hit : pixels_message->getData()) {
//     //     pixel_hits.push_back(&hit);
//     // }

//     // std::sort(pixel_hits.begin(), pixel_hits.end(), [](const PixelHit* a, const PixelHit* b) {
//     //     return a->getSignal() > b->getSignal();
//     // });
//     std::vector<const PixelHit*> sorted_pixel_hits(pixel_hits.begin(), pixel_hits.end());
//     std::sort(sorted_pixel_hits.begin(), sorted_pixel_hits.end(), [](const PixelHit* a, const PixelHit* b) {
//         return a->getSignal() > b->getSignal();
//     });

//     for(const auto& seed_candidate : sorted_pixel_hits) {
//         if(usedPixel[seed_candidate] || seed_candidate->getSignal() < seed_threshold_) {
//             continue;
//         }

//         Cluster cluster(seed_candidate);
//         usedPixel[seed_candidate] = true;
//         //LOG_STATUS.source("plot_AllpixData::doClustering") << "Creating new cluster with seed: " << seed_candidate->getPixel().getIndex();

//         std::vector<const PixelHit*> to_check;
//         to_check.push_back(seed_candidate);

//         while(!to_check.empty()) {
//             const PixelHit* current_pixel = to_check.front();
//             to_check.erase(to_check.begin());

//             const auto& current_index = current_pixel->getIndex();

//             for(const auto& neighbor_candidate : pixel_hits) {
//                 if(usedPixel[neighbor_candidate] || neighbor_candidate->getSignal() < neighbor_threshold_ || neighbor_candidate == current_pixel) {
//                     continue;
//                 }

//                 const auto& neighbor_index = neighbor_candidate->getIndex();

//                 if(std::abs(current_index.x() - neighbor_index.x()) <= 1 &&
//                    std::abs(current_index.y() - neighbor_index.y()) <= 1) {
//                     cluster.addPixelHit(neighbor_candidate);
//                     usedPixel[neighbor_candidate] = true;
//                     to_check.push_back(neighbor_candidate);
//                     //LOG_DEBUG.source("plot_AllpixData::doClustering") << "Adding pixel: " << neighbor_candidate->getPixel().getIndex();
//                 }
//             }

//             // for(const auto& neighbor_candidate : pixel_hits) {
//             //     if(usedPixel[neighbor_candidate] || neighbor_candidate->getSignal() < neighbor_threshold_ || neighbor_candidate == current_pixel) {
//             //         continue;
//             //     }

//                 // if(detector_->getModel()->areNeighbors(current_pixel->getIndex(), neighbor_candidate->getIndex(), 1)) {
//                 //     cluster.addPixelHit(neighbor_candidate);
//                 //     usedPixel[neighbor_candidate] =  true;
//                 //     to_check.push_back(neighbor_candidate);
//                 //     LOG_STATUS.source("plot_AllpixData::doClustering") << "Adding pixel: " << neighbor_candidate->getPixel().getIndex();
//                 // }
//             // }
//         }
//         clusters.push_back(cluster);

//     }
//     return clusters;
// }

// std::vector<const MCParticle*> plot_AllpixData::getPrimaryParticles(const std::vector<MCParticle*>& mc_particles) {
//     std::vector<const MCParticle*> primaries;

//     // // Loop over all MCParticles available
//     // for(const auto& mc_particle : mcparticle_message->getData()) {
//     //     // Check for possible parents:
//     //     const auto* parent = mc_particle.getParent();
//     //     if(parent != nullptr) {
//     //         LOG(TRACE) << "MCParticle " << mc_particle.getParticleID();
//     //         continue;
//     //     }

//     //     // This particle has no parent particles in the regarded sensor, return it.
//     //     LOG(TRACE) << "MCParticle " << mc_particle.getParticleID() << " (primary)";
//     //     primaries.push_back(&mc_particle);
//     // }

//     for(const auto& mc_particle : mc_particles) {
//         const auto* parent = mc_particle->getParent();
//         if(parent == nullptr) {
//             // This particle has no parent particles in the regarded sensor, return it.
//             LOG_STATUS.source("plot_AllpixData::getPrimaryParticles") << "MCParticle " << mc_particle->getParticleID() << " (primary)";
//             primaries.push_back(mc_particle);
//         }
//     }

//     return primaries;
// }



