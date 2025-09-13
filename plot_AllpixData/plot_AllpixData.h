// #ifndef PLOT_ALLPIXDATA_H
// #define PLOT_ALLPIXDATA_H

// #include <string>
// #include <vector>
// #include <utility>
// #include <functional>
// #include <cmath>
// #include <algorithm>
// #include <map>
// #include <iostream>

// // ROOT headers
// #include <TGraphErrors.h>
// #include <TMultiGraph.h>
// #include <TCanvas.h>
// #include <TLegend.h>
// #include <TLatex.h>
// #include <TFile.h>
// #include <TH1D.h>
// #include <TH2D.h>
// #include <TF1.h>
// #include <TStyle.h>
// #include <TColor.h>
// #include <TArrow.h>
// #include <TMarker.h>
// #include <TList.h>
// #include <TObjString.h>
// #include <TProfile2D.h>
// #include <TSystem.h>
// #include <TCollection.h>
// #include <TSystem.h>
// #include <TTree.h>

// // My function headers
// #include "Messenger/Messenger.h"
// #include "plot_histogram/plot_histogram.h"
// #include "plot_ExperimentData/plot_ExperimentData.h"

// // Other function headers
// #include "tools/json.hpp"
// #include "tools/cxxopts.hpp"
// #include "tools/Cluster.hpp"

// // Allpix headers
// //#include "/home/towa/package/allpix/allpix_dev/src/core/module/Module.hpp"
// #include "/home/towa/package/allpix/allpix_dev/src/objects/MCParticle.hpp"
// #include "/home/towa/package/allpix/allpix_dev/src/objects/PixelHit.hpp"
// #include "/home/towa/package/allpix/allpix_dev/src/objects/PropagatedCharge.hpp"
// //#include "/home/towa/package/allpix/allpix_dev/src/core/config/exceptions.h"


// enum class SimDataSource {
//     SingleChipSim,
//     TelescopeSim
// };

// namespace allpix {
//     class plot_AllpixData {
//     public:
//         plot_AllpixData();
//         ~plot_AllpixData();


//         void run_analysis();

//         std::vector<Cluster> doClustering(const std::vector<PixelHit*>& pixel_hits) const;

//         std::vector<const MCParticle*> getPrimaryParticles(const std::vector<MCParticle*>& mc_particles);

//     private:
//         //std::vector<PlotConfig> load_jsonConfigs(const nlohmann::json& j);
//         std::atomic<unsigned long> total_hits_{};
//         //std::shared_ptr<Detector> detector_;

//         ROOT::Math::XYVector matching_cut_{};

//         ROOT::Math::XYVector track_resolution_{};

//         //Histogram<TH2D> inPixel_hit_map, inPixel_mcHit_map;
//         //Histogram<TProfile2D> inPixel_residual, inPixel_cluster_charge, inPixel_seed_charge, inPixel_cluster_size, inPixel_electron_driftTime, inPixel_hole_driftTime;
        
//         double seed_threshold_;
//         double neighbor_threshold_;
//         //Histogram<TH1D> cluster_charge_histo_, cluster_size_histo_, seed_charge_histo_;

//         std::string file_name_{};
//         std::unique_ptr<TFile> data_file_{};
//         TTree* tree_{};
//         long long n_entries_{};

//         std::vector<PixelHit*>* pixel_hits_ = nullptr;
//         std::vector<MCParticle*>* mc_particles_ = nullptr;
//         std::vector<PropagatedCharge*>* propagated_charges_ = nullptr;
//     };
// }

// #endif // PLOT_ALLPIXDATA_H