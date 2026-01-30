import ROOT
# import uproot
import numpy as np
import math
from dataclasses import dataclass, field
from typing import List, Dict, Any
from ROOT import gSystem
from tqdm import tqdm
import argparse
import time
from datetime import datetime
import multiprocessing
import os
import glob
import traceback

gSystem.Load("/home/towa/package/allpix/install/lib/libAllpixObjects.so")
ROOT.ROOT.DisableImplicitMT()

# C++のPixelHit, MCParticleなどのデータ構造を模倣するためのデータクラス
@dataclass(slots=True)
class PixelHit:
    signal: float
    pixel_index_x: int
    pixel_index_y: int
    mc_particles: List[Any]

@dataclass(slots=True)
class MCParticle:
    particle_id: int
    parent: Any
    local_reference_point: np.ndarray

@dataclass(slots=True)
class PropagatedCharge:
    charge_type: str  # 'electron' or 'hole'
    global_time: float
    mc_particle: MCParticle
    local_creation_pos: np.ndarray

@dataclass(slots=True)
class Cluster:
    """Data class to hold cluster information"""
    seed_pixel_hit: PixelHit
    pixel_hits: List[PixelHit] = field(default_factory=list)
    
    def __post_init__(self):
        self.pixel_hits.append(self.seed_pixel_hit)

    def add_pixel_hit(self, hit: PixelHit):
        self.pixel_hits.append(hit)

    @property
    def charge(self) -> float:
        return sum(hit.signal for hit in self.pixel_hits)

    @property
    def size(self) -> int:
        return len(self.pixel_hits)
    
    def get_position(self, model, one_bit: bool) -> np.ndarray:
        pos_x, pos_y, pos_z = 0.0, 0.0, 0.0
        
        if one_bit:
            if self.size == 0:
                return model.get_pixel_center(
                    self.seed_pixel_hit.pixel_index_x,
                    self.seed_pixel_hit.pixel_index_y,
                )
            
            for hit in self.pixel_hits:
                center = model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
                pos_x += center[0]
                pos_y += center[1]
                pos_z += center[2]
            return np.array([pos_x / self.size, pos_y / self.size, pos_z / self.size])
        
        else:
            total_charge = self.charge
            if total_charge == 0:
                return model.get_pixel_center(
                    self.seed_pixel_hit.pixel_index_x,
                    self.seed_pixel_hit.pixel_index_y
                )

            for hit in self.pixel_hits:
                center = model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
                pos_x += center[0] * hit.signal
                pos_y += center[1] * hit.signal
                pos_z += center[2] * hit.signal
            return np.array([pos_x / total_charge, pos_y / total_charge, pos_z / total_charge])

    def get_mc_particles(self) -> List[MCParticle]:
        all_particles = []
        for hit in self.pixel_hits:
            all_particles.extend(hit.mc_particles)
        
        unique_particles = {}
        for particle in all_particles:
            if particle.particle_id not in unique_particles:
                unique_particles[particle.particle_id] = particle
        
        return list(unique_particles.values())


class DetectorModel:
    """Helper class to hold detector geometry information"""
    def __init__(self, pixel_size_x_um, pixel_size_y_um):
        self.pixel_size = np.array([pixel_size_x_um, pixel_size_y_um])

    def get_pixel_center(self, ix, iy):
        return np.array([(ix + 0.5) * self.pixel_size[0], (iy + 0.5) * self.pixel_size[1], 0.0])

    def are_neighbors(self, index1_x, index1_y, index2_x, index2_y, distance=1):
        return abs(index1_x - index2_x) <= distance and abs(index1_y - index2_y) <= distance

    def get_pixel_index(self, position):
        ix = math.floor(position[0] / self.pixel_size[0])
        iy = math.floor(position[1] / self.pixel_size[1])
        return ix, iy

class AnalysisPixelModule:
    """Python class reproducing C++ AnalysisPixelModule logic"""
    FOUR_NEIGHBORS = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    def __init__(self, config):
        self.config = config
        
        pitch_x = config.get("pixel_pitch_x", 22.5)
        pitch_y = config.get("pixel_pitch_y", 22.5)

        self.detector_name = config.get("detector_name", "detector")
        self.detector_model = DetectorModel(pixel_size_x_um=pitch_x, pixel_size_y_um=pitch_y)
        self.seed_threshold = config.get("seed_threshold", 1000)
        self.neighbor_threshold = config.get("neighbor_threshold", 500)
        self.histograms = {}

        self.one_bit_processing = config.get("one_bit", False)
        self.fill_3d = config.get("fill_3d", False)
        self.worker_id = config.get("worker_id", 0)
        self.max_cluster_size_hist = 10

        if self.one_bit_processing:
            self.neighbor_threshold = self.seed_threshold

        # Open file and TTrees using PyROOT
        self.input_file = ROOT.TFile.Open(config["file_name"])
        self.pixel_tree = self.input_file.Get("PixelHit")
        self.mcp_tree = self.input_file.Get("MCParticle")
        self.propagated_tree = self.input_file.Get("PropagatedCharge")

        if not self.pixel_tree or not self.mcp_tree or not self.propagated_tree:
            raise RuntimeError(f"File '{config['file_name']}' is missing required TTrees.")

        self.n_entries_total = self.pixel_tree.GetEntries()
        self.start_entry = config.get("start_entry", 0)
        self.end_entry = config.get("end_entry", self.n_entries_total)
        self.n_entries = self.end_entry - self.start_entry

        # Synchronize trees
        self.pixel_tree.AddFriend(self.mcp_tree)
        self.pixel_tree.AddFriend(self.propagated_tree)

        self.counter_names = {
            "Total Events": self.n_entries,
            "Skipped: No pixel hits": 0,
            "Skipped: Multiple primary particles": 0,
            "Clusters Checked": 0,
            "Skipped: Hit at EDGE events": 0,
            "Skipped: Cluster has no primary particle": 0,
            "Skipped: Residual > 40 um": 0,
            "Clusters Accepted": 0,
        }

        print(f"Worker processing {self.n_entries} events (range {self.start_entry} to {self.end_entry}) for file '{config['file_name']}'.")

    def setup_histograms(self):
        print("Creating histograms...")
        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        inpixel_bins_x = 50
        inpixel_bins_y = 50
        max_cluster_charge = self.config.get("max_cluster_charge", 20000)
        charge_bin = int(max_cluster_charge / 20)

        # Histogram definitions
        self.histograms["inPixel_neighbor_charge_sum"] = ROOT.TProfile2D("inPixel_neighbor_charge_sum", ";x/pitch [um];y/pitch [um];neighbor charge sum [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_efficiency"] = ROOT.TProfile2D("inPixel_efficiency", ";x/pitch [um];y/pitch [um];efficiency", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_cluster_size"] = ROOT.TProfile2D("inPixel_cluster_size", ";x/pitch [um];y/pitch [um];cluster size", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_cluster_charge"] = ROOT.TProfile2D("inPixel_cluster_charge", ";x/pitch [um];y/pitch [um];cluster charge [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_seed_charge"] = ROOT.TProfile2D("inPixel_seed_charge", ";x/pitch [um];y/pitch [um];seed charge [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_r"] = ROOT.TProfile2D("inPixel_residual_r", ";x/pitch [um];y/pitch [um];residual r [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_x"] = ROOT.TProfile2D("inPixel_residual_x", ";x/pitch [um];y/pitch [um];residual x [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_y"] = ROOT.TProfile2D("inPixel_residual_y", ";x/pitch [um];y/pitch [um];residual y [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_xy2"] = ROOT.TProfile2D("inPixel_residual_xy2", ";x/pitch [um];y/pitch [um];x + y / 2 [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["cluster_charge"] = ROOT.TH1D("cluster_charge", ";charge [ke];counts", charge_bin, 0, max_cluster_charge)
        self.histograms["cluster_size"] = ROOT.TH1D("cluster_size", ";cluster size;counts", 20, 0.5, 20.5)
        self.histograms["seed_charge"] = ROOT.TH1D("seed_charge", ";charge [ke];counts", charge_bin, 0, max_cluster_charge)
        self.histograms["residual_x"] = ROOT.TH1D("residual_x", ";residual x [um];counts", 10000, -40, 40)
        self.histograms["residual_y"] = ROOT.TH1D("residual_y", ";residual y [um];counts", 10000, -40, 40)
        self.histograms["residual_r"] = ROOT.TH1D("residual_r", ";residual r [um];counts", 10000, 0, 40)
        self.histograms["cluster_neighbor_charge_sum"] = ROOT.TH1D("cluster_neighbor_charge_sum", ";charge [ke];counts", charge_bin, 0, max_cluster_charge)
        self.histograms["cluster_neighbor_charge"] = ROOT.TH1D("cluster_neighbor_charge", ";charge [ke];counts", charge_bin, 0, max_cluster_charge)
        self.histograms["seedCharge_vs_clusterSize"] = ROOT.TH2D("seedCharge_vs_clusterSize", ";charge [ke];cluster size", charge_bin, 0, max_cluster_charge, 20, 0.5, 20.5)
        self.histograms["neighborChargeSum_vs_clusterSize"] = ROOT.TH2D("neighborChargeSum_vs_clusterSize", ";charge [ke];cluster size", charge_bin, 0, max_cluster_charge, 20, 0.5, 20.5)
        self.histograms["clusterCharge_vs_clusterSize"] = ROOT.TH2D("clusterCharge_vs_clusterSize", ";charge [ke];cluster size", charge_bin, 0, max_cluster_charge, 20, 0.5, 20.5)
        self.histograms["seedCharge_vs_neighborChargeSum"] = ROOT.TH2D("seedCharge_vs_neighborChargeSum", ";seed charge [ke];neighbor charge [ke]", charge_bin, 0, max_cluster_charge, 1000, 0, 20)
        self.histograms["clusterSize_vs_clusterCharge"] = ROOT.TH2D("clusterSize_vs_clusterCharge", "Cluster Size vs Cluster Charge;Cluster Charge [e];Cluster Size [pixels]", 100, 0, max_cluster_charge, 20, 0.5, 20.5)
        self.histograms["inPixel_seed_ratio"] = ROOT.TProfile2D("inPixel_seed_ratio", ";x/pitch [um];y/pitch [um];Seed Charge Ratio (Q_seed / Q_clus)", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_seed_charge_vs_neighborChargeSum"] = ROOT.TProfile2D("inPixel_seed_charge_vs_neighborChargeSum", ";x/pitch [um];y/pitch [um];seed charge / neighbor", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_seed_charge_vs_clusterSize"] = ROOT.TProfile2D("inPixel_seed_charge_vs_clusterSize", ";x/pitch [um];y/pitch [um];seed charge / cluster size", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_neighborChargeSum_vs_clusterSize"] = ROOT.TProfile2D("inPixel_neighborChargeSum_vs_clusterSize", ";x/pitch [um];y/pitch [um];neighbor / cluster size", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_cluster_charge_vs_clusterSize"] = ROOT.TProfile2D("inPixel_cluster_charge_vs_clusterSize", ";x/pitch [um];y/pitch [um];cluster_charge_vs_clusterSize", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_multi_hit_prob"] = ROOT.TProfile2D("inPixel_multi_hit_prob", ";x/pitch [um];y/pitch [um];Multi-Pixel Hit Probability", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)

        n_counters = len(self.counter_names)
        self.histograms["counters"] = ROOT.TH1D("counters", "Event Summary;category;counts", n_counters, 0, n_counters)
        for i, name in enumerate(self.counter_names):
            self.histograms["counters"].GetXaxis().SetBinLabel(i+1, name)

        for i in range(1, self.max_cluster_size_hist + 1):
            self.histograms[f"seed_charge_size_{i}"] = ROOT.TH1D(f"seed_charge_size_{i}", f"Cluster Seed Charge {i};charge [ke];counts", charge_bin, 0, max_cluster_charge)
            self.histograms[f"residual_x_size_{i}"] = ROOT.TH1D(f"residual_x_size_{i}", f"Residual X (Size {i});residual x [um];counts", 10000, -40, 40)
            self.histograms[f"residual_y_size_{i}"] = ROOT.TH1D(f"residual_y_size_{i}", f"Residual Y (Size {i});residual y [um];counts", 10000, -40, 40)
            self.histograms[f"residual_r_size_{i}"] = ROOT.TH1D(f"residual_r_size_{i}", f"Residual R (Size {i});residual r [um];counts", 10000, 0, 40)

        suffix_plus = f"{self.max_cluster_size_hist + 1}_plus"
        self.histograms[f"seed_charge_size_{suffix_plus}"] = ROOT.TH1D(f"seed_charge_size_{suffix_plus}", f"Seed Charge for Cluster Size > {self.max_cluster_size_hist};charge [ke];counts", charge_bin, 0, max_cluster_charge)
        self.histograms[f"residual_x_size_{suffix_plus}"] = ROOT.TH1D(f"residual_x_size_{suffix_plus}", f"Residual X (Size > {self.max_cluster_size_hist});residual x [um];counts", 10000, -40, 40)
        self.histograms[f"residual_y_size_{suffix_plus}"] = ROOT.TH1D(f"residual_y_size_{suffix_plus}", f"Residual Y (Size > {self.max_cluster_size_hist});residual y [um];counts", 10000, -40, 40)
        self.histograms[f"residual_r_size_{suffix_plus}"] = ROOT.TH1D(f"residual_r_size_{suffix_plus}", f"Residual R (Size > {self.max_cluster_size_hist});residual r [um];counts", 10000, 0, 40)

        self.histograms["inPixel_electron_driftTime_90p"] = ROOT.TProfile2D("inPixel_electron_driftTime_90p", ";x/pitch [um];y/pitch [um];90% electron drift time [ns]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["electron_driftTime_90p"] = ROOT.TH1D("electron_driftTime_90p", ";90% electron drift time [ns];counts", 20000, 0, 2)
        self.histograms["drift_time_vs_distance"] = ROOT.TH2D("drift_time_vs_distance", ";distance to pixel center [um];drift time [ns];counts", 200, 0, 100, 1000, 0, 0.2)
        self.histograms["cluster_charge_vs_drift_time"] = ROOT.TH2D("cluster_charge_vs_drift_time", ";Cluster Charge [ke];90% electron drift time [ns];Counts", charge_bin, 0, max_cluster_charge, 1000, 0, 0.2)
        self.histograms["cluster_size_vs_drift_time"] = ROOT.TH2D("cluster_size_vs_drift_time", ";Cluster Size [pixels];90% electron drift time [ns];Counts", 20, 0.5, 20.5, 1000, 0, 0.2)
        self.histograms["drift_time_spectrum"] = ROOT.TH1D("drift_time_spectrum", "All Electron Drift Time;Drift Time [ns];Counts", 20000, 0, 2)
        self.histograms["drift_time_vs_depth"] = ROOT.TH2D("drift_time_vs_depth", "Drift Time vs Creation Depth;Creation Z [um];drift time [ns]", 100, -30, 30, 20000, 0, 100)

        if self.fill_3d:
            half_pitch_x = pitch_x / 2.0
            half_pitch_y = pitch_y / 2.0
            sensor_half_thickness = 25.0
            self.histograms["drift_time_map_xyz"] = ROOT.TProfile3D("drift_time_map_xyz", ";In-Pixel X [um];In-Pixel Y [um];Depth Z [um];Average Drift Time [ns]", 10, -half_pitch_x, half_pitch_x, 10, -half_pitch_y, half_pitch_y, 25, -sensor_half_thickness, sensor_half_thickness)
        
        self.histograms["drift_time_seed"] = ROOT.TH1D("drift_time_seed", ";drift time (seed pixel) [ns];counts (per carrier)", 1000, 0, 0.2)
        self.histograms["drift_time_neighbor"] = ROOT.TH1D("drift_time_neighbor", ";drift time (neighbor pixel) [ns];counts (per carrier)", 1000, 0, 0.2)

    def do_clustering(self, pixel_hits: List[PixelHit]) -> List[Cluster]:
        clusters = []
        used_pixels = set()
        hit_map = {(p.pixel_index_x, p.pixel_index_y): p for p in pixel_hits}
        sorted_pixel_hits = sorted(pixel_hits, key=lambda p: p.signal, reverse=True)

        for seed_candidate in sorted_pixel_hits:
            seed_id = (seed_candidate.pixel_index_x, seed_candidate.pixel_index_y)
            if seed_id in used_pixels or seed_candidate.signal < self.seed_threshold:
                continue

            cluster = Cluster(seed_candidate)
            used_pixels.add(seed_id)
            to_check_queue = [seed_candidate]

            while to_check_queue:
                current_pixel = to_check_queue.pop(0)
                for dx, dy in self.FOUR_NEIGHBORS:
                    neighbor_id = (current_pixel.pixel_index_x + dx, current_pixel.pixel_index_y + dy)
                    if neighbor_id in hit_map and neighbor_id not in used_pixels:
                        neighbor_hit = hit_map[neighbor_id]
                        if neighbor_hit.signal >= self.neighbor_threshold:
                            cluster.add_pixel_hit(neighbor_hit)
                            used_pixels.add(neighbor_id)
                            to_check_queue.append(neighbor_hit)

            clusters.append(cluster)
        return clusters

    def get_primary_particles(self, mc_particles: Dict[int, MCParticle]) -> List[MCParticle]:
        return [p for p in mc_particles.values() if p.parent is None]

    def run_analysis(self):
        """Execute event loop and fill histograms directly (no batching)"""
        self.setup_histograms()
        branch_name = self.config.get("detector_name", "CE65")

        N_PIXELS_X, N_PIXELS_Y = 48, 24
        MARGIN = 2

        self.counters = {name: 0 for name in self.counter_names}
        self.counters["Total Events"] = self.n_entries

        pitch_x, pitch_y = self.detector_model.pixel_size[0], self.detector_model.pixel_size[1]
        offset_x, offset_y = pitch_x / 2, pitch_y / 2
        
        # Cache histogram objects for performance
        h_eff = self.histograms["inPixel_efficiency"]
        h_clus_charge = self.histograms["cluster_charge"]
        h_clus_size = self.histograms["cluster_size"]
        h_seed_charge = self.histograms["seed_charge"]
        h_res_x = self.histograms["residual_x"]
        h_res_y = self.histograms["residual_y"]
        h_res_r = self.histograms["residual_r"]
        h_ip_clus_size = self.histograms["inPixel_cluster_size"]
        h_ip_res_r = self.histograms["inPixel_residual_r"]
        h_ip_res_x = self.histograms["inPixel_residual_x"]
        h_ip_res_y = self.histograms["inPixel_residual_y"]
        h_ip_res_xy2 = self.histograms["inPixel_residual_xy2"]
        h_ip_seed_charge = self.histograms["inPixel_seed_charge"]
        h_ip_clus_charge = self.histograms["inPixel_cluster_charge"]
        h_ip_neighbor_charge_sum = self.histograms["inPixel_neighbor_charge_sum"]
        h_size_vs_charge = self.histograms["clusterSize_vs_clusterCharge"]
        h_ip_seed_ratio = self.histograms["inPixel_seed_ratio"]
        h_ip_multi_hit = self.histograms["inPixel_multi_hit_prob"]
        h_ip_seed_vs_neighbor = self.histograms["inPixel_seed_charge_vs_neighborChargeSum"]
        h_ip_seed_vs_size = self.histograms["inPixel_seed_charge_vs_clusterSize"]
        h_ip_neighbor_vs_size = self.histograms["inPixel_neighborChargeSum_vs_clusterSize"]
        h_ip_clus_vs_size = self.histograms["inPixel_cluster_charge_vs_clusterSize"]
        h_neighbor_sum = self.histograms["cluster_neighbor_charge_sum"]
        h_neighbor_sum_vs_size = self.histograms["neighborChargeSum_vs_clusterSize"]
        h_seed_vs_neighbor = self.histograms["seedCharge_vs_neighborChargeSum"]
        h_seed_vs_size = self.histograms["seedCharge_vs_clusterSize"]
        h_clus_vs_size = self.histograms["clusterCharge_vs_clusterSize"]
        h_neighbor_charge = self.histograms["cluster_neighbor_charge"]
        
        h_drift_time_vs_dist = self.histograms.get("drift_time_vs_distance")
        h_drift_time_map_xyz = self.histograms.get("drift_time_map_xyz") if self.fill_3d else None
        h_drift_time_seed = self.histograms.get("drift_time_seed")
        h_drift_time_neighbor = self.histograms.get("drift_time_neighbor")
        h_drift_time_spectrum = self.histograms.get("drift_time_spectrum")
        h_drift_time_vs_depth = self.histograms.get("drift_time_vs_depth")
        h_dt_90p = self.histograms["electron_driftTime_90p"]
        h_ip_dt_90p = self.histograms["inPixel_electron_driftTime_90p"]
        h_clus_charge_vs_dt = self.histograms["cluster_charge_vs_drift_time"]
        h_clus_size_vs_dt = self.histograms["cluster_size_vs_drift_time"]

        tqdm_iterator = tqdm(
            range(self.start_entry, self.end_entry),
            desc=f"Worker {self.worker_id} (Events {self.start_entry}-{self.end_entry})", 
            position=self.worker_id,
            leave=False
        )

        for i_global in tqdm_iterator:
            self.pixel_tree.GetEntry(i_global)

            # --- 1. MCParticle Data Decoding ---
            mcp_objects = getattr(self.mcp_tree, branch_name)
            mc_particles: Dict[int, MCParticle] = {}
            for mcp_obj in mcp_objects:
                pid = mcp_obj.GetUniqueID()
                parent_ptr = mcp_obj.getParent()
                parent_id = parent_ptr.GetUniqueID() if parent_ptr else 0
                ref_point = mcp_obj.getLocalReferencePoint()
                mc_particles[pid] = MCParticle(
                    particle_id=pid, parent=parent_id,
                    local_reference_point=np.array([ref_point.X()*1000, ref_point.Y()*1000, ref_point.Z()*1000])
                )

            for mcp in mc_particles.values():
                parent_id = mcp.parent
                mcp.parent = mc_particles.get(parent_id) if parent_id != 0 else None

            # --- 2. Primary Particle and ROI Check ---
            primary_particles = self.get_primary_particles(mc_particles)
            if len(primary_particles) != 1:
                self.counters["Skipped: Multiple primary particles"] += 1
                continue
            
            target_particle = primary_particles[0]
            particle_pos = target_particle.local_reference_point + np.array([offset_x, offset_y, 0])
            ix_mc, iy_mc = self.detector_model.get_pixel_index(particle_pos)
            is_in_roi = (MARGIN <= ix_mc < N_PIXELS_X - MARGIN and MARGIN <= iy_mc < N_PIXELS_Y - MARGIN)
            pixel_center_mc = self.detector_model.get_pixel_center(ix_mc, iy_mc)
            in_pixel_pos_mc = particle_pos - pixel_center_mc

            # --- 3. PixelHit processing and Clustering ---
            pixel_objects = getattr(self.pixel_tree, branch_name)
            is_reconstructed_successfully = False
            best_cluster, residual_vec, residual_r = None, None, None

            if len(pixel_objects) > 0:
                pixel_hits = []
                for hit_obj in pixel_objects:
                    hit_mc_particles = [mc_particles[p.GetUniqueID()] for p in hit_obj.getMCParticles() if p.GetUniqueID() in mc_particles]
                    pixel_hits.append(PixelHit(
                        signal=hit_obj.getSignal(),
                        pixel_index_x=hit_obj.getPixel().getIndex().X(),
                        pixel_index_y=hit_obj.getPixel().getIndex().Y(),
                        mc_particles=hit_mc_particles
                    ))
                clusters = self.do_clustering(pixel_hits)
                if clusters:
                    best_cluster = max(clusters, key=lambda c: c.seed_pixel_hit.signal)
                    self.counters["Clusters Checked"] += 1
                    cluster_particles_ids = {p.particle_id for p in best_cluster.get_mc_particles()}
                    if target_particle.particle_id in cluster_particles_ids:
                        cluster_pos = best_cluster.get_position(self.detector_model, one_bit=self.one_bit_processing)
                        residual_vec = particle_pos - cluster_pos
                        residual_r = np.linalg.norm(residual_vec[:2])
                        if residual_r <= 40:
                            is_reconstructed_successfully = True

            # --- 4. Efficiency Fill (Directly) ---
            if is_in_roi:
                h_eff.Fill(in_pixel_pos_mc[0], in_pixel_pos_mc[1], 1.0 if is_reconstructed_successfully else 0.0)

            if not is_reconstructed_successfully:
                if len(pixel_objects) == 0: self.counters["Skipped: No pixel hits"] += 1
                elif best_cluster is None: pass
                elif residual_vec is None: self.counters["Skipped: Cluster has no primary particle"] += 1
                else: self.counters["Skipped: Residual > 40 um"] += 1
                continue

            # --- 5. Analysis of Accepted Clusters ---
            seed_x, seed_y = best_cluster.seed_pixel_hit.pixel_index_x, best_cluster.seed_pixel_hit.pixel_index_y
            if not (MARGIN <= seed_x < N_PIXELS_X - MARGIN and MARGIN <= seed_y < N_PIXELS_Y - MARGIN):
                self.counters["Skipped: Hit at EDGE events"] += 1
                continue

            self.counters["Clusters Accepted"] += 1
            clus = best_cluster
            
            # Drift Time Prep
            propagated_objects = getattr(self.propagated_tree, branch_name, [])
            charge_map: Dict[int, List[PropagatedCharge]] = {}
            for prop_obj in propagated_objects:
                mc_particle_ptr = prop_obj.getMCParticle()
                if not mc_particle_ptr: continue
                parent_mcp = mc_particles.get(mc_particle_ptr.GetUniqueID())
                if not parent_mcp: continue
                
                carrier_type = prop_obj.getType()
                charge_type = 'electron' if carrier_type == 255 else 'hole' if carrier_type == 1 else None
                if not charge_type: continue

                pos = prop_obj.getLocalPosition()
                charge_map.setdefault(parent_mcp.particle_id, []).append(PropagatedCharge(
                    charge_type=charge_type, global_time=prop_obj.getGlobalTime(),
                    mc_particle=parent_mcp, local_creation_pos=np.array([pos.X()*1000, pos.Y()*1000, pos.Z()*1000])
                ))

            # In-Pixel Position Re-calculation
            ix, iy = self.detector_model.get_pixel_index(particle_pos)
            pixel_center = self.detector_model.get_pixel_center(ix, iy)
            in_pixel_pos = particle_pos - pixel_center
            
            non_seed_hits = [hit for hit in clus.pixel_hits if hit is not clus.seed_pixel_hit]
            sum_non_seed_charge = sum(hit.signal for hit in non_seed_hits)

            # Direct Histogram Fills
            h_clus_charge.Fill(clus.charge)
            h_clus_size.Fill(clus.size)
            h_seed_charge.Fill(clus.seed_pixel_hit.signal)
            h_res_x.Fill(residual_vec[0])
            h_res_y.Fill(residual_vec[1])
            h_res_r.Fill(residual_r)
            h_ip_clus_size.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.size)
            h_ip_res_r.Fill(in_pixel_pos[0], in_pixel_pos[1], residual_r)
            h_ip_res_x.Fill(in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[0]))
            h_ip_res_y.Fill(in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[1]))
            h_ip_res_xy2.Fill(in_pixel_pos[0], in_pixel_pos[1], (abs(residual_vec[0]) + abs(residual_vec[1])) / 2)
            h_ip_seed_charge.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal)
            h_ip_clus_charge.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.charge)
            h_ip_neighbor_charge_sum.Fill(in_pixel_pos[0], in_pixel_pos[1], sum_non_seed_charge)
            h_size_vs_charge.Fill(clus.charge, clus.size)
            
            if clus.charge > 0:
                h_ip_seed_ratio.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal / clus.charge)

            h_ip_multi_hit.Fill(in_pixel_pos[0], in_pixel_pos[1], 1.0 if clus.size > 1 else 0.0)

            if sum_non_seed_charge > 0:
                h_ip_seed_vs_neighbor.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal / sum_non_seed_charge)
                h_neighbor_sum.Fill(sum_non_seed_charge)
                h_neighbor_sum_vs_size.Fill(sum_non_seed_charge, clus.size)
                h_seed_vs_neighbor.Fill(clus.seed_pixel_hit.signal, sum_non_seed_charge)
            
            h_ip_seed_vs_size.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal / clus.size)
            h_ip_neighbor_vs_size.Fill(in_pixel_pos[0], in_pixel_pos[1], sum_non_seed_charge / clus.size)
            h_ip_clus_vs_size.Fill(in_pixel_pos[0], in_pixel_pos[1], clus.charge / clus.size)
            h_seed_vs_size.Fill(clus.seed_pixel_hit.signal, clus.size)
            h_clus_vs_size.Fill(clus.charge, clus.size)

            for hit in non_seed_hits:
                h_neighbor_charge.Fill(hit.signal)

            # Per Size Histograms
            size_suffix = f"_size_{clus.size}" if clus.size <= self.max_cluster_size_hist else f"_size_{self.max_cluster_size_hist+1}_plus"
            self.histograms[f"seed_charge{size_suffix}"].Fill(clus.seed_pixel_hit.signal)
            self.histograms[f"residual_x{size_suffix}"].Fill(residual_vec[0])
            self.histograms[f"residual_y{size_suffix}"].Fill(residual_vec[1])
            self.histograms[f"residual_r{size_suffix}"].Fill(residual_r)

            # Drift Time
            electron_drift_times_ps = []
            for hit in clus.pixel_hits:
                pixel_center_pos = self.detector_model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
                for mcp in hit.mc_particles:
                    for pc in charge_map.get(mcp.particle_id, []):
                        if pc.charge_type == 'electron':
                            dt_ns = pc.global_time * 1e-3
                            electron_drift_times_ps.append(pc.global_time)
                            dist_um = np.linalg.norm(pc.local_creation_pos - pixel_center_pos)
                            
                            if h_drift_time_vs_dist: h_drift_time_vs_dist.Fill(dist_um, dt_ns)
                            if h_drift_time_spectrum: h_drift_time_spectrum.Fill(dt_ns)
                            if h_drift_time_vs_depth: h_drift_time_vs_depth.Fill(pc.local_creation_pos[2], dt_ns)
                            if h_drift_time_map_xyz:
                                rel = pc.local_creation_pos - pixel_center_pos
                                h_drift_time_map_xyz.Fill(rel[0], rel[1], rel[2], dt_ns)
                            
                            if hit is clus.seed_pixel_hit:
                                if h_drift_time_seed: h_drift_time_seed.Fill(dt_ns)
                            else:
                                if h_drift_time_neighbor: h_drift_time_neighbor.Fill(dt_ns)

            if electron_drift_times_ps:
                electron_drift_times_ps.sort()
                time_90p_ns = electron_drift_times_ps[int(len(electron_drift_times_ps) * 0.9)] * 1e-3
                h_dt_90p.Fill(time_90p_ns)
                h_ip_dt_90p.Fill(in_pixel_pos[0], in_pixel_pos[1], time_90p_ns)
                h_clus_charge_vs_dt.Fill(clus.charge, time_90p_ns)
                h_clus_size_vs_dt.Fill(clus.size, time_90p_ns)

        print("Analysis finished.")
        print("----------------- Report -----------------")
        for key, value in self.counters.items():
            print(f"{key:<40}:{value}")
        print("------------------------------------------")

    def finalize(self):
        for i, name in enumerate(self.counter_names):
            count = self.counters.get(name, 0)
            if self.histograms.get("counters"):
                self.histograms["counters"].SetBinContent(i + 1, count)

        output_filename = self.config.get("output_file_name")
        output_file = ROOT.TFile(output_filename, "RECREATE")

        dir_inpixel = output_file.mkdir("inPixel")
        dir_per_size = output_file.mkdir("PerClusterSize")
        dir_drifttime = output_file.mkdir("DriftTime")

        dir_size_specific = {}
        for i in range(1, self.max_cluster_size_hist + 1):
            dir_size_specific[f"_size_{i}"] = dir_per_size.mkdir(f"clsize_{i}")
        suffix_plus = f"{self.max_cluster_size_hist + 1}_plus"
        dir_size_specific[f"_size_{suffix_plus}"] = dir_per_size.mkdir(f"clsize_{suffix_plus}")

        for name, histo in self.histograms.items():
            written = False
            for suffix, directory in dir_size_specific.items():
                if name.endswith(suffix):
                    directory.cd()
                    histo.Write()
                    written = True
                    break
            if written: continue

            if name.startswith("inPixel_"):
                dir_inpixel.cd()
                histo.Write()
            elif (name.startswith("drift_time_") or name.startswith("electron_driftTime_") or 
                  name.endswith("_vs_drift_time") or name.endswith("_vs_depth")):
                dir_drifttime.cd()
                histo.Write()
            else:
                output_file.cd()
                histo.Write()

        output_file.Close()
        print(f"Histograms written to {output_filename}")

def run_worker(args):
    start_index, end_index, worker_id, config = args
    config["output_file_name"] = f"analysis_py_part_{worker_id}.root"
    config["start_entry"], config["end_entry"], config["worker_id"] = start_index, end_index, worker_id
    try:
        analyzer = AnalysisPixelModule(config)
        analyzer.run_analysis()
        analyzer.finalize()
        return config["output_file_name"]
    except Exception as e:
        print(f"--- Worker {worker_id} FAILED ---")
        traceback.print_exc()
        return None

if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output", default="analysis_py.root")
    parser.add_argument("-st", "--seed_threshold", type=int, default=0)
    parser.add_argument("-nt", "--neighbor_threshold", type=int, default=0)
    parser.add_argument("-p", "--pixel_pitch", type=float, default=22.5)
    parser.add_argument("-v", "--voltage", type=int, default=10)
    parser.add_argument("-b", "--beam_type", default="e3GeV")
    parser.add_argument("-m", "--model", default="masetti")
    parser.add_argument("-n", "--name", default="CE65")
    parser.add_argument("--one_bit", action="store_true")
    parser.add_argument("--fill3D", action="store_true")
    parser.add_argument("-j", "--cores", type=int, default=1)
    args = parser.parse_args()

    config = {
        "file_name": args.input_file, "output_file_name": args.output,
        "detector_name": args.name, "pixel_pitch_x": args.pixel_pitch, "pixel_pitch_y": args.pixel_pitch,
        "seed_threshold": args.seed_threshold, "neighbor_threshold": args.neighbor_threshold,
        "max_cluster_charge": 20000, "one_bit": args.one_bit, "fill_3d": args.fill3D
    }

    temp_file = ROOT.TFile.Open(config["file_name"])
    total_entries = temp_file.Get("PixelHit").GetEntries()
    temp_file.Close()

    chunk_size = math.ceil(total_entries / args.cores)
    tasks = [(i*chunk_size, min((i+1)*chunk_size, total_entries), i, config.copy()) for i in range(args.cores) if i*chunk_size < total_entries]

    with multiprocessing.Pool(processes=args.cores) as pool:
        results = pool.map(run_worker, tasks)
        partial_files = [res for res in results if res]

    if partial_files:
        os.system(f"hadd -f {config['output_file_name']} {' '.join(partial_files)}")
        for f in partial_files: os.remove(f)

    elapsed = time.time() - start_time
    print(f"Total time: {elapsed:.2f} seconds.")