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
import gc

gSystem.Load("/home/towa/package/allpix/install/lib/libAllpixObjects.so")
ROOT.ROOT.DisableImplicitMT()

# C++のPixelHit, MCParticleなどのデータ構造を模倣するためのデータクラス
# (この部分は変更ありません)
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
    """クラスタ情報を保持するデータクラス"""
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
        # pos_x, pos_y, pos_z = 0.0, 0.0, 0.0
        
        # if one_bit:
        #     if self.size == 0:
        #         return model.get_pixel_center(
        #             self.seed_pixel_hit.pixel_index_x,
        #             self.seed_pixel_hit.pixel_index_y,
        #         )
            
        #     for hit in self.pixel_hits:
        #         center = model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
        #         pos_x += center[0]
        #         pos_y += center[1]
        #         pos_z += center[2]
        #     return np.array([pos_x / self.size, pos_y / self.size, pos_z / self.size])
        
        # else:
        #     total_charge = self.charge
        #     if total_charge == 0:
        #         return model.get_pixel_center(
        #             self.seed_pixel_hit.pixel_index_x,
        #             self.seed_pixel_hit.pixel_index_y
        #         )

        #     for hit in self.pixel_hits:
        #         center = model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
        #         pos_x += center[0] * hit.signal
        #         pos_y += center[1] * hit.signal
        #         pos_z += center[2] * hit.signal
        #     return np.array([pos_x / total_charge, pos_y / total_charge, pos_z / total_charge])
        """Calculate cluster position using numpy vector operations"""
        # Collect all hit coordinates into an (N, 3) array
        centers = np.array([model.get_pixel_center(h.pixel_index_x, h.pixel_index_y) for h in self.pixel_hits])
        
        if one_bit:
            # Average of pixel centers
            return np.mean(centers, axis=0)
        else:
            # Charge-weighted average
            signals = np.array([h.signal for h in self.pixel_hits])
            total_charge = np.sum(signals)
            if total_charge <= 0:
                return centers[0]
            # Use np.average for weighted mean: sum(centers * signals) / sum(signals)
            return np.average(centers, axis=0, weights=signals)

    def get_mc_particles(self) -> List[MCParticle]:
        # all_particles = []
        # for hit in self.pixel_hits:
        #     all_particles.extend(hit.mc_particles)
        
        # # particle_idを使って重複を削除し、順番を保持する
        # unique_particles = {}
        # for particle in all_particles:
        #     if particle.particle_id not in unique_particles:
        #         unique_particles[particle.particle_id] = particle
        
        # return list(unique_particles.values())
        unique_particles = {p.particle_id: p for hit in self.pixel_hits for p in hit.mc_particles}
        return list(unique_particles.values())


class DetectorModel:
    """検出器のジオメトリ情報を保持するヘルパークラス"""
    def __init__(self, pixel_size_x_um, pixel_size_y_um):
        self.pixel_size = np.array([pixel_size_x_um, pixel_size_y_um])

    def get_pixel_center(self, ix, iy):
        return np.array([(ix + 0.5) * self.pixel_size[0], (iy + 0.5) * self.pixel_size[1], 0.0])

    def are_neighbors(self, index1_x, index1_y, index2_x, index2_y, distance=1):
        return abs(index1_x - index2_x) <= distance and abs(index1_y - index2_y) <= distance

    def get_pixel_index(self, position_um: np.ndarray):
        #ix = int(position[0] / self.pixel_size[0])
        #iy = int(position[1] / self.pixel_size[1])
        # #return ix, iy
        # ix = math.floor(position[0] / self.pixel_size[0])
        # iy = math.floor(position[1] / self.pixel_size[1])
        # return ix, iy
        """Convert position to pixel indices using numpy floor"""
        indices = np.floor(position_um[:2] / self.pixel_size)
        return int(indices[0]), int(indices[1])

class AnalysisPixelModule:
    FOUR_NEIGHBORS = np.array([(0, 1), (0, -1), (1, 0), (-1, 0)])

    def __init__(self, config):
        # self.config = config
        
        # pitch_x = config.get("pixel_pitch_x", 22.5)
        # pitch_y = config.get("pixel_pitch_y", 22.5)

        # self.detector_name = config.get("detector_name", "detector")
        # self.detector_model = DetectorModel(pixel_size_x_um=pitch_x, pixel_size_y_um=pitch_y)
        # self.seed_threshold = config.get("seed_threshold", 1000)
        # self.neighbor_threshold = config.get("neighbor_threshold", 500)
        # self.histograms = {}

        # self.one_bit_processing = config.get("one_bit", False)

        # self.fill_3d = config.get("fill_3d", False)

        # self.worker_id = config.get("worker_id", 0)

        # self.max_cluster_size_hist = 10

        # if self.one_bit_processing:
        #     self.neighbor_threshold = self.seed_threshold

        # # PyROOTでファイルとTTreeを開く
        # self.input_file = ROOT.TFile.Open(config["file_name"])
        # self.pixel_tree = self.input_file.Get("PixelHit")
        # self.mcp_tree = self.input_file.Get("MCParticle")
        # self.n_entries = self.pixel_tree.GetEntries()
        # #self.n_entries = 100  
        # self.propagated_tree = self.input_file.Get("PropagatedCharge")

        # if not self.pixel_tree or not self.mcp_tree or not self.propagated_tree:
        #     raise RuntimeError(f"File '{config['file_name']}' is missing required TTrees (PixelHit, MCParticle, or PropagatedCharge).")

        # self.n_entries_total = self.pixel_tree.GetEntries()
        # self.start_entry = config.get("start_entry", 0)
        # # configでend_entryが指定されなければ、全イベントを処理
        # self.end_entry = config.get("end_entry", self.n_entries_total)
        # # このワーカーが処理するイベント数
        # self.n_entries = self.end_entry - self.start_entry

        # # 2つのTTreeが同じイベントを指すように同期させる
        # self.pixel_tree.AddFriend(self.mcp_tree)
        # self.pixel_tree.AddFriend(self.propagated_tree)

        # self.counter_names = {
        #     "Total Events": self.n_entries,
        #     "Skipped: No pixel hits": 0,
        #     "Skipped: Multiple primary particles": 0,
        #     "Clusters Checked": 0,
        #     "Skipped: Hit at EDGE events": 0,
        #     "Skipped: Cluster has no primary particle": 0,
        #     "Skipped: Residual > 40 um": 0,
        #     "Clusters Accepted": 0,
        # }
        self.config = config
        pitch = np.array([config.get("pixel_pitch_x", 22.5), config.get("pixel_pitch_y", 22.5)])
        self.detector_model = DetectorModel(pitch[0], pitch[1])
        self.seed_threshold = config.get("seed_threshold", 1000)
        self.neighbor_threshold = config.get("neighbor_threshold", 500)
        self.histograms = {}
        self.one_bit_processing = config.get("one_bit", False)
        self.fill_3d = config.get("fill_3d", False)
        self.worker_id = config.get("worker_id", 0)
        if self.one_bit_processing:
            self.neighbor_threshold = self.seed_threshold

        self.input_file = ROOT.TFile.Open(config["file_name"])
        self.pixel_tree = self.input_file.Get("PixelHit")
        self.mcp_tree = self.input_file.Get("MCParticle")
        self.propagated_tree = self.input_file.Get("PropagatedCharge")
        self.pixel_tree.AddFriend(self.mcp_tree)
        self.pixel_tree.AddFriend(self.propagated_tree)

        self.start_entry = config.get("start_entry", 0)
        self.end_entry = config.get("end_entry", self.pixel_tree.GetEntries())
        self.max_cluster_size_hist = 10
        self.counter_names = [
            "Total Events",
            "Skipped: No pixel hits",
            "Skipped: Multiple primary particles",
            "Clusters Checked",
            "Skipped: Hit at EDGE events",
            "Skipped: Cluster has no primary particle",
            "Skipped: Residual > 40 um",
            "Clusters Accepted"
        ]
        self.counters = {name: 0 for name in self.counter_names}
        # self.counters = {name: 0 for name in [
        #     "Total Events",
        #     "Skipped: No pixel hits",
        #     "Skipped: Multiple primary particles",
        #     "Clusters Checked",
        #     "Skipped: Hit at EDGE events",
        #     "Skipped: Cluster has no primary particle",
        #     "Skipped: Residual > 40 um",
        #     "Clusters Accepted"
        # ]}

        #self.four_neighbors = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        self.n_entries = self.end_entry - self.start_entry
        print(f"Worker processing {self.n_entries} events (range {self.start_entry} to {self.end_entry}) for file '{config['file_name']}'.")
        print(f"File '{config['file_name']}' opened with {self.n_entries} events.")

    def setup_histograms(self):
        # (この関数は変更ありません)
        print("Creating histograms...")
        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        inpixel_bins_x = 50
        inpixel_bins_y = 50
        max_cluster_charge = self.config.get("max_cluster_charge", 20000)
        charge_bin = int(max_cluster_charge / 20)
        self.histograms["inPixel_neighbor_charge_sum"] = ROOT.TProfile2D("inPixel_neighbor_charge_sum", ";x/pitch [um];y/pitch [um];neighbor charge sum [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_efficiency"] = ROOT.TProfile2D("inPixel_efficiency", ";x/pitch [um];y/pitch [um];efficiency", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_cluster_size"] = ROOT.TProfile2D("inPixel_cluster_size", ";x/pitch [um];y/pitch [um];cluster size", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_cluster_charge"] = ROOT.TProfile2D("inPixel_cluster_charge", ";x/pitch [um];y/pitch [um];cluster charge [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_seed_charge"] = ROOT.TProfile2D("inPixel_seed_charge", ";x/pitch [um];y/pitch [um];seed charge [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_r"] = ROOT.TProfile2D("inPixel_residual_r", ";x/pitch [um];y/pitch [um];residual r [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_x"] = ROOT.TProfile2D("inPixel_residual_x", ";x/pitch [um];y/pitch [um];residual x [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_y"] = ROOT.TProfile2D("inPixel_residual_y", ";x/pitch [um];y/pitch [um];residual y [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_xy2"] = ROOT.TProfile2D("inPixel_residual_xy2", ";x/pitch [um];y/pitch [um];x + y / 2 [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual"] = ROOT.TProfile2D("inPixel_residual", ";x/pitch [um];y/pitch [um];residual [um]", inpixel_bins_x, -1000, 1000, inpixel_bins_y, -1000, 1000)
        #self.histograms["test"] = ROOT.TH2D("test", "test", 1000, -50, 50, 1000, -50, 50)
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
        
        self.histograms["clusterSize_vs_clusterCharge"] = ROOT.TH2D(
            "clusterSize_vs_clusterCharge",
            "Cluster Size vs Cluster Charge;Cluster Charge [e];Cluster Size [pixels]",
            100, 0, max_cluster_charge,  # X軸: 電荷 (0 ~ 20000e)
            20, 0.5, 20.5                # Y軸: サイズ (1 ~ 20)
        )

        # self.histograms["prof_clusterSize_vs_clusterCharge"] = ROOT.TProfile2D(
        #     "prof_clusterSize_vs_clusterCharge",
        #     "Average Cluster Size vs Cluster Charge;Cluster Charge [e];Average Cluster Size [pixels]",
        #     100, 0, max_cluster_charge,
        #     20, 0.5, 20.5
        # )

        self.histograms["inPixel_seed_ratio"] = ROOT.TProfile2D(
            "inPixel_seed_ratio", 
            ";x/pitch [um];y/pitch [um];Seed Charge Ratio (Q_seed / Q_clus)", 
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )

        self.histograms["inPixel_seed_charge_vs_neighborChargeSum"] = ROOT.TProfile2D(
            "inPixel_seed_charge_vs_neighborChargeSum",
            ";x/pitch [um];y/pitch [um];seed charge / neighbor",
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )

        self.histograms["inPixel_seed_charge_vs_clusterSize"] = ROOT.TProfile2D(
            "inPixel_seed_charge_vs_clusterSize",
            ";x/pitch [um];y/pitch [um];seed charge / cluster size",
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )

        self.histograms["inPixel_neighborChargeSum_vs_clusterSize"] = ROOT.TProfile2D(
            "inPixel_neighborChargeSum_vs_clusterSize",
            ";x/pitch [um];y/pitch [um];neighbor / cluster size",
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )

        self.histograms["inPixel_cluster_charge_vs_clusterSize"] = ROOT.TProfile2D(
            "inPixel_cluster_charge_vs_clusterSize",
            ";x/pitch [um];y/pitch [um];cluster_charge_vs_clusterSize",
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )
        
        self.histograms["inPixel_multi_hit_prob"] = ROOT.TProfile2D(
            "inPixel_multi_hit_prob", 
            ";x/pitch [um];y/pitch [um];Multi-Pixel Hit Probability", 
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )

        self.percentiles = self.config.get("percentiles", [50, 60, 70, 80, 90])
        for p in self.percentiles:
            h_name = f"electron_driftTime_{p}p"
            self.histograms[h_name] = ROOT.TH1D(h_name, f";{p}% electron drift time [ns];counts", 2000, 0, 2)
            prof_name = f"inPixel_electron_driftTime_{p}p"
            self.histograms[prof_name] = ROOT.TProfile2D(
                prof_name,
                f";x/pitch [um];y/pitch [um];{p}% electron drift time [ns]",
                inpixel_bins_x, -pitch_x / 2, pitch_x / 2,
                inpixel_bins_y, -pitch_y / 2, pitch_y / 2
            )
            self.histograms[f"inPixel_driftTime_seed_{p}p"] = ROOT.TProfile2D(
                f"inPixel_driftTime_seed_{p}p", 
                ";x/pitch [um];y/pitch [um];average drift time (seed) [ns]", 
                inpixel_bins_x, -pitch_x/2, pitch_x/2, inpixel_bins_y, -pitch_y/2, pitch_y/2
            )
            self.histograms[f"inPixel_driftTime_neighbor_{p}p"] = ROOT.TProfile2D(
                f"inPixel_driftTime_neighbor_{p}p", 
                ";x/pitch [um];y/pitch [um];average drift time (neighbor) [ns]", 
                inpixel_bins_x, -pitch_x/2, pitch_x/2, inpixel_bins_y, -pitch_y/2, pitch_y/2
            )

        self.target_carrier_counts = self.config.get("target_carrier_counts", [48, 96, 144, 150, 160, 170, 140, 180, 190,200, 240, 312])
        for c in self.target_carrier_counts:
            h_name_tc = f"driftTime_to_{c}e"
            self.histograms[h_name_tc] = ROOT.TH1D(h_name_tc, f";time to reach {c}e [ns];counts", 2000, 0, 2)

            prof_name_tc = f"inPixel_driftTime_to_{c}e"
            self.histograms[prof_name_tc] = ROOT.TProfile2D(
                prof_name_tc,
                f";x/pitch [um];y/pitch [um];time to reach {c}e [ns]",
                inpixel_bins_x, -pitch_x/2, pitch_x/2,
                inpixel_bins_y, -pitch_y/2, pitch_y/2
            )


        # self.histograms["inPixel_driftTime_seed_carriers"] = ROOT.TProfile2D(
        #     "inPixel_driftTime_seed_carriers", 
        #     ";x/pitch [um];y/pitch [um];average drift time (seed) [ns]", 
        #     inpixel_bins_x, -pitch_x/2, pitch_x/2, inpixel_bins_y, -pitch_y/2, pitch_y/2
        # )
        # self.histograms["inPixel_driftTime_neighbor_carriers"] = ROOT.TProfile2D(
        #     "inPixel_driftTime_neighbor_carriers", 
        #     ";x/pitch [um];y/pitch [um];average drift time (neighbor) [ns]", 
        #     inpixel_bins_x, -pitch_x/2, pitch_x/2, inpixel_bins_y, -pitch_y/2, pitch_y/2
        # )

        n_counters = len(self.counter_names)
        self.histograms["counters"] = ROOT.TH1D("counters", "Event Summary;category;counts", n_counters, 0, n_counters)
        for i, name in enumerate(self.counter_names):
            self.histograms["counters"].GetXaxis().SetBinLabel(i+1, name)

        for i in range(1, self.max_cluster_size_hist + 1):
            hist_name_seed = f"seed_charge_size_{i}"
            hist_title_seed = f"Cluster Seed Charge {i};charge [ke];counts"
            self.histograms[hist_name_seed] = ROOT.TH1D(hist_name_seed, hist_title_seed, charge_bin, 0, max_cluster_charge)
            
            self.histograms[f"residual_x_size_{i}"] = ROOT.TH1D(f"residual_x_size_{i}", f"Residual X (Size {i});residual x [um];counts", 10000, -40, 40)
            self.histograms[f"residual_y_size_{i}"] = ROOT.TH1D(f"residual_y_size_{i}", f"Residual Y (Size {i});residual y [um];counts", 10000, -40, 40)
            self.histograms[f"residual_r_size_{i}"] = ROOT.TH1D(f"residual_r_size_{i}", f"Residual R (Size {i});residual r [um];counts", 10000, 0, 40)

        suffix_plus = f"{self.max_cluster_size_hist + 1}_plus"
        hist_name_seed_large = f"seed_charge_size_{suffix_plus}"
        self.histograms[hist_name_seed_large] = ROOT.TH1D(hist_name_seed_large, f"Seed Charge for Cluster Size > {self.max_cluster_size_hist};charge [ke];counts", charge_bin, 0, max_cluster_charge)

        # NEW: Residuals overflow
        self.histograms[f"residual_x_size_{suffix_plus}"] = ROOT.TH1D(f"residual_x_size_{suffix_plus}", f"Residual X (Size > {self.max_cluster_size_hist});residual x [um];counts", 10000, -40, 40)
        self.histograms[f"residual_y_size_{suffix_plus}"] = ROOT.TH1D(f"residual_y_size_{suffix_plus}", f"Residual Y (Size > {self.max_cluster_size_hist});residual y [um];counts", 10000, -40, 40)
        self.histograms[f"residual_r_size_{suffix_plus}"] = ROOT.TH1D(f"residual_r_size_{suffix_plus}", f"Residual R (Size > {self.max_cluster_size_hist});residual r [um];counts", 10000, 0, 40)

        #self.histograms["inPixel_electron_driftTime_90p"] = ROOT.TProfile2D("inPixel_electron_driftTime_90p", ";x/pitch [um];y/pitch [um];90% electron drift time [ns]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        #self.histograms["electron_driftTime_90p"] = ROOT.TH1D("electron_driftTime_90p", ";90% electron drift time [ns];counts", 20000, 0, 2)

        self.histograms["drift_time_vs_distance"] = ROOT.TH2D(
            "drift_time_vs_distance",
            ";distance to pixel center [um];drift time [ns];counts",
            200, 0, 100,
            1000, 0, 0.2
        )
        self.histograms["cluster_charge_vs_drift_time"] = ROOT.TH2D(
            "cluster_charge_vs_drift_time",
            ";Cluster Charge [ke];90% electron drift time [ns];Counts",
            charge_bin, 0, max_cluster_charge,
            1000, 0, 0.2
        )
        self.histograms["cluster_size_vs_drift_time"] = ROOT.TH2D(
            "cluster_size_vs_drift_time",
            ";Cluster Size [pixels];90% electron drift time [ns];Counts",
            20, 0.5, 20.5,
            1000, 0, 0.2
        )

        self.histograms["drift_time_spectrum"] = ROOT.TH1D("drift_time_spectrum", "All Electron Drift Time;Drift Time [ns];Counts", 20000, 0, 2)
        self.histograms["drift_time_vs_depth"] = ROOT.TH2D("drift_time_vs_depth", "Drift Time vs Creation Depth;Creation Z [um];drift time [ns]", 100, -30, 30, 20000, 0, 100)

        if self.fill_3d:
            half_pitch_x = self.detector_model.pixel_size[0] / 2.0
            half_pitch_y = self.detector_model.pixel_size[1] / 2.0
            sensor_half_thickness = 25.0 # センサ厚 50umを想定

            self.histograms["drift_time_map_xyz"] = ROOT.TProfile3D(
                "drift_time_map_xyz",
                ";In-Pixel X [um];In-Pixel Y [um];Depth Z [um];Average Drift Time [ns]",
                10, -half_pitch_x, half_pitch_x,
                10, -half_pitch_y, half_pitch_y,
                25, -sensor_half_thickness, sensor_half_thickness
            )
        
        self.histograms["drift_time_seed"] = ROOT.TH1D(
            "drift_time_seed",
            ";drift time (seed pixel) [ns];counts (per carrier)",
            1000, 0, 0.2
        )
        self.histograms["drift_time_neighbor"] = ROOT.TH1D(
            "drift_time_neighbor",
            ";drift time (neighbor pixel) [ns];counts (per carrier)",
            1000, 0, 0.2
        )

    def do_clustering(self, pixel_hits: List[PixelHit]) -> List[Cluster]:
        clusters = []
        used_pixels = set()
        hit_map = {(p.pixel_index_x, p.pixel_index_y): p for p in pixel_hits}
        sorted_hits = sorted(pixel_hits, key=lambda p: p.signal, reverse=True)

        for seed in sorted_hits:
            seed_id = (seed.pixel_index_x, seed.pixel_index_y)
            if seed_id in used_pixels or seed.signal < self.seed_threshold:
                continue

            cluster = Cluster(seed)
            used_pixels.add(seed_id)
            queue = [seed]

            while queue:
                curr = queue.pop(0)
                for dx, dy in self.FOUR_NEIGH_OFFSET:
                    nb_id = (curr.pixel_index_x + dx, curr.pixel_index_y + dy)
                    if nb_id in hit_map and nb_id not in used_pixels:
                        nb_hit = hit_map[nb_id]
                        if nb_hit.signal >= self.neighbor_threshold:
                            cluster.add_pixel_hit(nb_hit)
                            used_pixels.add(nb_id)
                            queue.append(nb_hit)
            clusters.append(cluster)
        return clusters

    @property
    def FOUR_NEIGH_OFFSET(self):
        return [(0, 1), (0, -1), (1, 0), (-1, 0)]

    def get_primary_particles(self, mc_particles: Dict[int, MCParticle]) -> List[MCParticle]:
        # (この関数は変更ありません)
        return [p for p in mc_particles.values() if p.parent is None]

    def _fill_histograms_from_buffer(self, buffer: Dict[str, list]):
        """バッファに溜まったデータからヒストグラムを更新する"""
        for key, histo in self.histograms.items():
            data_list = buffer.get(key)
            if not data_list: continue
            
            if isinstance(histo, ROOT.TProfile3D):
                for x, y, z, val in data_list:
                    histo.Fill(x, y, z, val)
            elif isinstance(histo, ROOT.TProfile2D):
                for x, y, val in data_list:
                    histo.Fill(x, y, val)
            elif isinstance(histo, ROOT.TH2D):
                for x_val, y_val in data_list:
                    histo.Fill(x_val, y_val)
            elif isinstance(histo, ROOT.TH1D):
                for val in data_list:
                    histo.Fill(val)
            data_list.clear()
        

    # def run_analysis(self):
    #     """PyROOTを使ってイベントループを実行し、解析を行う"""
    #     self.setup_histograms()
        
    #     branch_name = "CE65"

    #     N_PIXELS_X = 48
    #     N_PIXELS_Y = 24
    #     MARGIN = 2

    #     self.counters = {name: 0 for name in self.counter_names}
    #     self.counters["Total Events"] = self.n_entries

    #     pitch_x = self.detector_model.pixel_size[0]
    #     pitch_y = self.detector_model.pixel_size[1]
    #     offset_x = pitch_x / 2
    #     offset_y = pitch_y / 2
        
    #     BATCH_SIZE = 100
    #     buffer = {hist_name: [] for hist_name in self.histograms}

    #     # 高速化のために変数をキャッシュ
    #     h_drift_time_vs_dist = self.histograms.get("drift_time_vs_distance")
    #     h_drift_time_map_xyz = self.histograms.get("drift_time_map_xyz") if self.fill_3d else None
    #     h_drift_time_seed = self.histograms.get("drift_time_seed")
    #     h_drift_time_neighbor = self.histograms.get("drift_time_neighbor")
    #     h_drift_time_spectrum = self.histograms.get("drift_time_spectrum")
    #     h_drift_time_vs_depth = self.histograms.get("drift_time_vs_depth")

    #     tqdm_iterator = tqdm(
    #         range(self.start_entry, self.end_entry),
    #         desc=f"Worker {self.worker_id} (Events {self.start_entry}-{self.end_entry})", 
    #         position=self.worker_id,
    #         leave=False
    #     )

    #     for i_global in tqdm_iterator:
    #         self.pixel_tree.GetEntry(i_global)

    #         # --- 1. MCParticle データデコード (位置を最初に特定するため移動) ---
    #         mcp_objects = getattr(self.mcp_tree, branch_name)
    #         mc_particles: Dict[int, MCParticle] = {}
            
    #         for mcp_obj in mcp_objects:
    #             pid = mcp_obj.GetUniqueID()
    #             parent_ptr = mcp_obj.getParent()
    #             parent_id = parent_ptr.GetUniqueID() if parent_ptr else 0
    #             ref_point = mcp_obj.getLocalReferencePoint()
    #             mc_particles[pid] = MCParticle(
    #                 particle_id=pid,
    #                 parent=parent_id,
    #                 local_reference_point=np.array([ref_point.X()*1000, ref_point.Y()*1000, ref_point.Z()*1000])
    #             )

    #         for mcp in mc_particles.values():
    #             parent_id = mcp.parent
    #             if parent_id != 0 and parent_id in mc_particles:
    #                 mcp.parent = mc_particles[parent_id]
    #             else:
    #                 mcp.parent = None

    #         # --- 2. Primary Particle チェックとROI判定 ---
    #         primary_particles = self.get_primary_particles(mc_particles)
            
    #         if len(primary_particles) != 1:
    #             self.counters["Skipped: Multiple primary particles"] += 1
    #             continue
            
    #         # Primary粒子の位置を特定
    #         target_particle = primary_particles[0]
    #         particle_pos = target_particle.local_reference_point + np.array([offset_x, offset_y, 0])
            
    #         # MC粒子がROI内にあるか確認 (Efficiencyの分母)
    #         ix_mc, iy_mc = self.detector_model.get_pixel_index(particle_pos)
            
    #         is_in_roi = (MARGIN <= ix_mc < N_PIXELS_X - MARGIN and MARGIN <= iy_mc < N_PIXELS_Y - MARGIN)
            
    #         pixel_center_mc = self.detector_model.get_pixel_center(ix_mc, iy_mc)
    #         in_pixel_pos_mc = particle_pos - pixel_center_mc

    #         # --- 3. PixelHit 処理とクラスタリング ---
    #         pixel_objects = getattr(self.pixel_tree, branch_name)
            
    #         # Efficiency計算用のフラグ
    #         is_reconstructed_successfully = False
    #         best_cluster = None
    #         residual_vec = None
    #         residual_r = None

    #         # ヒットがある場合のみクラスタリング処理を行う
    #         if len(pixel_objects) > 0:
    #             pixel_hits = []
    #             for hit_obj in pixel_objects:
    #                 hit_mc_particles = []
    #                 for linked_mcp_ptr in hit_obj.getMCParticles():
    #                     linked_pid = linked_mcp_ptr.GetUniqueID()
    #                     if linked_pid in mc_particles:
    #                         hit_mc_particles.append(mc_particles[linked_pid])
                    
    #                 pixel_hits.append(PixelHit(
    #                     signal=hit_obj.getSignal(),
    #                     pixel_index_x=hit_obj.getPixel().getIndex().X(),
    #                     pixel_index_y=hit_obj.getPixel().getIndex().Y(),
    #                     mc_particles=hit_mc_particles
    #                 ))

    #             clusters = self.do_clustering(pixel_hits)

    #             if clusters:
    #                 best_cluster = max(clusters, key=lambda c: c.seed_pixel_hit.signal)
    #                 self.counters["Clusters Checked"] += 1
                    
    #                 # クラスタがPrimary粒子と紐付いているか確認
    #                 cluster_particles_ids = {p.particle_id for p in best_cluster.get_mc_particles()}
    #                 if target_particle.particle_id in cluster_particles_ids:
                        
    #                     # 残差計算
    #                     cluster_pos = best_cluster.get_position(self.detector_model, one_bit=self.one_bit_processing)
    #                     residual_vec = particle_pos - cluster_pos
    #                     residual_r = np.linalg.norm(residual_vec[:2])

    #                     cut_um = 40
    #                     if residual_r <= cut_um:
    #                         is_reconstructed_successfully = True

    #         # --- 4. Efficiency Plot の Fill ---
    #         # ROI内に入射したイベントであれば、検出できた(1)か否(0)かを記録
    #         if is_in_roi:
    #             eff_val = 1.0 if is_reconstructed_successfully else 0.0
    #             buffer["inPixel_efficiency"].append((in_pixel_pos_mc[0], in_pixel_pos_mc[1], eff_val))

    #         # --- 5. 検出失敗ならここで終了 (Efficiency 0 のケース) ---
    #         if not is_reconstructed_successfully:
    #             if len(pixel_objects) == 0:
    #                 self.counters["Skipped: No pixel hits"] += 1
    #             elif best_cluster is None:
    #                 # ヒットはあるがクラスタにならなかった（閾値など）
    #                 pass 
    #             elif residual_vec is None:
    #                 # クラスタはあるがPrimaryと紐付かなかった
    #                 self.counters["Skipped: Cluster has no primary particle"] += 1
    #             else:
    #                 # 残差カット落ち
    #                 self.counters["Skipped: Residual > 40 um"] += 1
    #             continue

    #         # --- 6. 以下、検出成功 (Efficiency 1) イベントの詳細解析 ---
    #         # ここから先は `best_cluster` が有効で、かつResidualカットも通過していることが保証される
            
    #         # Seed PixelがEDGEにある場合は除く（従来の解析ロジックを踏襲）
    #         # 注意: EfficiencyMapはMC位置基準でROIを切ったが、詳細プロットはSeed位置基準で切るのが通例
    #         seed_x = best_cluster.seed_pixel_hit.pixel_index_x
    #         seed_y = best_cluster.seed_pixel_hit.pixel_index_y

    #         if not (MARGIN <= seed_x < N_PIXELS_X - MARGIN and MARGIN <= seed_y < N_PIXELS_Y - MARGIN):
    #             self.counters["Skipped: Hit at EDGE events"] += 1
    #             continue

    #         self.counters["Clusters Accepted"] += 1
    #         clus = best_cluster # 変数名を合わせる

    #         # Propagated Charge の処理 (Drift Time用)
    #         propagated_objects = getattr(self.propagated_tree, branch_name, [])
    #         charge_map: Dict[int, List[PropagatedCharge]] = {}
    #         for prop_obj in propagated_objects:
    #             mc_particle_ptr = prop_obj.getMCParticle()
    #             if not mc_particle_ptr: continue
    #             parent_mcp = mc_particles.get(mc_particle_ptr.GetUniqueID())
    #             if not parent_mcp: continue
                
    #             carrier_type_enum = prop_obj.getType()
    #             if carrier_type_enum == 255: charge_type = 'electron'
    #             elif carrier_type_enum == 1: charge_type = 'hole'
    #             else: continue

    #             pos = prop_obj.getLocalPosition()
    #             creation_pos_um = np.array([pos.X() * 1000, pos.Y() * 1000, pos.Z() * 1000])
                
    #             prop_charge = PropagatedCharge(
    #                 charge_type=charge_type,
    #                 global_time=prop_obj.getGlobalTime(),
    #                 mc_particle=parent_mcp,
    #                 local_creation_pos=creation_pos_um
    #             )
    #             charge_map.setdefault(parent_mcp.particle_id, []).append(prop_charge)


    #         # In-Pixel Positionの再計算 (Seed基準: プロットの一貫性のため)
    #         ix, iy = self.detector_model.get_pixel_index(particle_pos)
    #         pixel_center = self.detector_model.get_pixel_center(ix, iy)
    #         in_pixel_pos = particle_pos - pixel_center
            
    #         non_seed_hits = [hit for hit in clus.pixel_hits if hit is not clus.seed_pixel_hit]
    #         sum_non_seed_charge = sum(hit.signal for hit in non_seed_hits)

    #         # ヒストグラムへのFill
    #         buffer["cluster_charge"].append(clus.charge)
    #         buffer["cluster_size"].append(clus.size)
    #         buffer["seed_charge"].append(clus.seed_pixel_hit.signal)
    #         buffer["residual_x"].append(residual_vec[0])
    #         buffer["residual_y"].append(residual_vec[1])
    #         buffer["residual_r"].append(residual_r)
    #         buffer["inPixel_cluster_size"].append((in_pixel_pos[0], in_pixel_pos[1], clus.size))
    #         buffer["inPixel_residual_r"].append((in_pixel_pos[0], in_pixel_pos[1], residual_r))
    #         buffer["inPixel_residual_x"].append((in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[0])))
    #         buffer["inPixel_residual_y"].append((in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[1])))
    #         buffer["inPixel_residual_xy2"].append((in_pixel_pos[0], in_pixel_pos[1], (abs(residual_vec[0]) + abs(residual_vec[1])) / 2))
    #         buffer["inPixel_seed_charge"].append((in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal))
    #         buffer["inPixel_cluster_charge"].append((in_pixel_pos[0], in_pixel_pos[1], clus.charge))
    #         buffer["inPixel_neighbor_charge_sum"].append((in_pixel_pos[0], in_pixel_pos[1], sum_non_seed_charge))
    #         buffer["clusterSize_vs_clusterCharge"].append((clus.charge, clus.size))
    #         #buffer["prof_clusterSize_vs_clusterCharge"].append((clus.charge, clus.size))
           
    #         # 1. Seed Charge Ratio (0.0 ~ 1.0)
    #         if clus.charge > 0:
    #             ratio = clus.seed_pixel_hit.signal / clus.charge
    #             buffer["inPixel_seed_ratio"].append((in_pixel_pos[0], in_pixel_pos[1], ratio))

    #         # 2. Multi-Pixel Hit Probability (0 or 1)
    #         # サイズが1より大きければ 1.0 (True), そうでなければ 0.0 (False)
    #         is_multi = 1.0 if clus.size > 1 else 0.0
    #         buffer["inPixel_multi_hit_prob"].append((in_pixel_pos[0], in_pixel_pos[1], is_multi))

    #         # inPixel_seed_charge_vs_neighborChargeSum
    #         if sum_non_seed_charge > 0:
    #             val_seed_vs_neighbor = clus.seed_pixel_hit.signal / sum_non_seed_charge
    #             buffer["inPixel_seed_charge_vs_neighborChargeSum"].append((in_pixel_pos[0], in_pixel_pos[1], val_seed_vs_neighbor))
            
    #         # inPixel_seed_charge_vs_clusterSize
    #         val_seed_vs_size = clus.seed_pixel_hit.signal / clus.size
    #         buffer["inPixel_seed_charge_vs_clusterSize"].append((in_pixel_pos[0], in_pixel_pos[1], val_seed_vs_size))

    #         # inPixel_neighborChargeSum_vs_clusterSize
    #         val_neighbor_vs_size = sum_non_seed_charge / clus.size
    #         buffer["inPixel_neighborChargeSum_vs_clusterSize"].append((in_pixel_pos[0], in_pixel_pos[1], val_neighbor_vs_size))

    #         # inPixel_cluster_charge_vs_clusterSize
    #         val_cluster_vs_size = clus.charge / clus.size
    #         buffer["inPixel_cluster_charge_vs_clusterSize"].append((in_pixel_pos[0], in_pixel_pos[1], val_cluster_vs_size))

    #         if sum_non_seed_charge > 0:
    #             buffer["cluster_neighbor_charge_sum"].append(sum_non_seed_charge)
    #             buffer["neighborChargeSum_vs_clusterSize"].append((sum_non_seed_charge, clus.size))
    #             buffer["seedCharge_vs_neighborChargeSum"].append((clus.seed_pixel_hit.signal, sum_non_seed_charge))

    #         buffer["seedCharge_vs_clusterSize"].append((clus.seed_pixel_hit.signal, clus.size))
    #         buffer["clusterCharge_vs_clusterSize"].append((clus.charge, clus.size))

    #         for hit in non_seed_hits:
    #             buffer["cluster_neighbor_charge"].append(hit.signal)

    #         size = clus.size
    #         seed_charge_ke = clus.seed_pixel_hit.signal

    #         if size <= self.max_cluster_size_hist:
    #             buffer[f"seed_charge_size_{size}"].append(seed_charge_ke)
    #             buffer[f"residual_x_size_{size}"].append(residual_vec[0])
    #             buffer[f"residual_y_size_{size}"].append(residual_vec[1])
    #             buffer[f"residual_r_size_{size}"].append(residual_r)
    #         else:
    #             suffix_plus = f"{self.max_cluster_size_hist + 1}_plus"
    #             buffer[f"seed_charge_size_{suffix_plus}"].append(seed_charge_ke)
    #             buffer[f"residual_x_size_{suffix_plus}"].append(residual_vec[0])
    #             buffer[f"residual_y_size_{suffix_plus}"].append(residual_vec[1])
    #             buffer[f"residual_r_size_{suffix_plus}"].append(residual_r)
            
    #         # --- Drift Time Analysis ---
    #         electron_drift_times_ps = []
            
    #         for hit in clus.pixel_hits:
    #             pixel_center_pos = self.detector_model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
                
    #             for mcp in hit.mc_particles:
    #                 propagated_charges = charge_map.get(mcp.particle_id, [])
                    
    #                 for pc in propagated_charges:
    #                     if pc.charge_type == 'electron':
    #                         drift_time_ns = pc.global_time * 1e-3
    #                         electron_drift_times_ps.append(pc.global_time)
                            
    #                         creation_pos = pc.local_creation_pos
    #                         distance_um = np.linalg.norm(creation_pos - pixel_center_pos)
                            
    #                         if h_drift_time_vs_dist: h_drift_time_vs_dist.Fill(distance_um, drift_time_ns)
    #                         if h_drift_time_spectrum: h_drift_time_spectrum.Fill(drift_time_ns)
    #                         if h_drift_time_vs_depth: h_drift_time_vs_depth.Fill(creation_pos[2], drift_time_ns)

    #                         if h_drift_time_map_xyz:
    #                             relative_pos = creation_pos - pixel_center_pos
    #                             h_drift_time_map_xyz.Fill(relative_pos[0], relative_pos[1], relative_pos[2], drift_time_ns)
                            
    #                         if hit is clus.seed_pixel_hit:
    #                             if h_drift_time_seed: h_drift_time_seed.Fill(drift_time_ns)
    #                         else:
    #                             if h_drift_time_neighbor: h_drift_time_neighbor.Fill(drift_time_ns)

    #         if electron_drift_times_ps:
    #             electron_drift_times_ps.sort()
    #             idx_90 = int(len(electron_drift_times_ps) * 0.9)
    #             time_90p_ns = electron_drift_times_ps[idx_90] * 1e-3

    #             buffer["electron_driftTime_90p"].append(time_90p_ns)
    #             buffer["inPixel_electron_driftTime_90p"].append((in_pixel_pos[0], in_pixel_pos[1], time_90p_ns))
    #             buffer["cluster_charge_vs_drift_time"].append((clus.charge, time_90p_ns))
    #             buffer["cluster_size_vs_drift_time"].append((clus.size, time_90p_ns))

    #         if (i_global + 1) % BATCH_SIZE == 0:
    #             self._fill_histograms_from_buffer(buffer)
        
    #     if any(buffer.values()):
    #         self._fill_histograms_from_buffer(buffer)

    #     print("Analysis finished.")
    #     print("----------------- Report -----------------")
    #     for key, value in self.counters.items():
    #         print(f"{key:<40}:{value}")
    #     print("------------------------------------------")
        
    # --- 変更箇所: 解析ループ内のメモリ管理を徹底 ---

    def run_analysis(self, events_per_chunk=500): # チャンクをさらに小さく設定可能に
        self.setup_histograms()
        branch = self.config.get("detector_name", "CE65")
        N_X, N_Y, MARGIN = 48, 24, 2
        
        BATCH_SIZE = 100
        buffer = {k: [] for k in self.histograms}
        
        pitch = self.detector_model.pixel_size
        offset_vec = np.array([pitch[0]/2, pitch[1]/2, 0.0])

        # 高速化のためにメソッドをローカル変数にキャッシュ
        get_pixel_center = self.detector_model.get_pixel_center
        get_pixel_index = self.detector_model.get_pixel_index

        for i_ev in tqdm(range(self.start_entry, self.end_entry), position=self.worker_id, leave=False):
            # 1. データの読み込み
            if self.pixel_tree.GetEntry(i_ev) <= 0: continue
            self.counters["Total Events"] += 1

            # --- 1. MCParticles (最小限の生成) ---
            mcp_objs = getattr(self.mcp_tree, branch)
            mc_particles = {}
            primaries = []
            for obj in mcp_objs:
                pid = obj.GetUniqueID()
                ref = obj.getLocalReferencePoint()
                # 必要な情報だけ保持
                mcp = MCParticle(pid, obj.getParent().GetUniqueID() if obj.getParent() else 0, 
                               np.array([ref.X(), ref.Y(), ref.Z()]) * 1000.0)
                mc_particles[pid] = mcp
                if mcp.parent == 0: primaries.append(mcp)

            if len(primaries) != 1:
                self.counters["Skipped: Multiple primary particles"] += 1
                # メモリ解放
                del mc_particles, primaries
                continue
            
            target_mcp = primaries[0]
            particle_pos = target_mcp.local_reference_point + offset_vec
            ix_mc, iy_mc = get_pixel_index(particle_pos)
            
            is_in_roi = (MARGIN <= ix_mc < N_X - MARGIN and MARGIN <= iy_mc < N_Y - MARGIN)
            in_pixel_pos_mc = particle_pos - get_pixel_center(ix_mc, iy_mc)

            # --- 2. PixelHits & Clustering ---
            pixel_objs = getattr(self.pixel_tree, branch)
            hits = []
            for o in pixel_objs:
                linked_mcp = [mc_particles[p.GetUniqueID()] for p in o.getMCParticles() if p.GetUniqueID() in mc_particles]
                hits.append(PixelHit(o.getSignal(), o.getPixel().getIndex().X(), o.getPixel().getIndex().Y(), linked_mcp))

            clusters = self.do_clustering(hits)
            best_cl = max(clusters, key=lambda c: c.seed_pixel_hit.signal) if clusters else None
            
            is_rec = False
            res_vec = None
            if best_cl:
                if target_mcp.particle_id in {p.particle_id for p in best_cl.get_mc_particles()}:
                    res_vec = particle_pos - best_cl.get_position(self.detector_model, self.one_bit_processing)
                    if np.linalg.norm(res_vec[:2]) <= 40: is_rec = True

            if is_in_roi:
                buffer["inPixel_efficiency"].append((in_pixel_pos_mc[0], in_pixel_pos_mc[1], 1.0 if is_rec else 0.0))

            # 検出失敗時の早期クリーンアップ
            if not is_rec:
                if not hits: self.counters["Skipped: No pixel hits"] += 1
                elif best_cl is None: pass
                elif res_vec is None: self.counters["Skipped: Cluster has no primary particle"] += 1
                else: self.counters["Skipped: Residual > 40 um"] += 1
                # 明示的な削除
                del hits, mc_particles, primaries, clusters
                continue

            # --- 3. Accepted Cluster Detail Analysis ---
            seed = best_cl.seed_pixel_hit
            if not (MARGIN <= seed.pixel_index_x < N_X - MARGIN and MARGIN <= seed.pixel_index_y < N_Y - MARGIN):
                self.counters["Skipped: Hit at EDGE events"] += 1
                del hits, mc_particles, primaries, clusters
                continue

            self.counters["Clusters Accepted"] += 1
            pix_center = get_pixel_center(seed.pixel_index_x, seed.pixel_index_y)
            in_pix = particle_pos - pix_center
            q_sum, q_seed, cl_size = best_cl.charge, seed.signal, best_cl.size
            res_r = np.linalg.norm(res_vec[:2])

            # Buffer Fill (省略せず記述が必要ですがメモリ負荷は低い)
            # Fill Buffer
            q_nb_sum = q_sum - q_seed
            in_x, in_y = in_pix[0], in_pix[1]
            buffer["cluster_charge"].append(q_sum)
            buffer["cluster_size"].append(cl_size)
            buffer["seed_charge"].append(q_seed)
            buffer["residual_x"].append(res_vec[0])
            buffer["residual_y"].append(res_vec[1])
            buffer["residual_r"].append(res_r)
            buffer["inPixel_cluster_size"].append((in_x, in_y, cl_size))
            buffer["inPixel_seed_charge"].append((in_x, in_y, q_seed))
            buffer["inPixel_cluster_charge"].append((in_x, in_y, q_sum))
            buffer["inPixel_neighbor_charge_sum"].append((in_x, in_y, q_nb_sum))
            buffer["inPixel_residual_r"].append((in_x, in_y, res_r))

            if q_sum > 0:
                buffer["inPixel_seed_ratio"].append((in_x, in_y, q_seed / q_sum))
            buffer["inPixel_multi_hit_prob"].append((in_x, in_y, 1.0 if cl_size > 1 else 0.0))

            # Cluster Size Specific
            s_suffix = f"_size_{cl_size}" if cl_size <= self.max_cluster_size_hist else f"_size_{self.max_cluster_size_hist + 1}_plus"
            buffer[f"seed_charge{s_suffix}"].append(q_seed)
            buffer[f"residual_x{s_suffix}"].append(res_vec[0])
            buffer[f"residual_y{s_suffix}"].append(res_vec[1])
            buffer[f"residual_r{s_suffix}"].append(res_r)

            # --- 4. Propagated Charge (ここが最重要: データクラスを作らない) ---
            # 直接ROOTオブジェクトから値を読み取ってFillすることでPythonのオーバーヘッドを無くす
            prop_objs = getattr(self.propagated_tree, branch)
            dts = []
            seed_mcp_ids = {p.particle_id for p in seed.mc_particles}

            all_dts, seed_dts, neighbor_dts = [], [], []

            seed_x = seed.pixel_index_x
            seed_y = seed.pixel_index_y

            for p_obj in prop_objs:
                if p_obj.getType() != 255: continue # Electrons
                mcp_ptr = p_obj.getMCParticle()
                if not mcp_ptr or mcp_ptr.GetUniqueID() != target_mcp.particle_id: continue
                
                dt = p_obj.getGlobalTime() # ps
                dt_ns = dt * 1e-3
                all_dts.append(dt)
                self.histograms["drift_time_spectrum"].Fill(dt_ns)

                pos = p_obj.getLocalPosition()
                cur_pix_x, cur_pix_y = get_pixel_index(np.array([pos.X()*1000 + offset_vec[0], 
                                                                pos.Y()*1000 + offset_vec[1]]))

                # if mcp_ptr.GetUniqueID() in seed_mcp_ids:
                #     seed_dts.append(dt)
                #     self.histograms["drift_time_seed"].Fill(dt_ns)
                # else:
                #     neighbor_dts.append(dt)
                #     self.histograms["drift_time_neighbor"].Fill(dt_ns)
                if cur_pix_x == seed_x and cur_pix_y == seed_y:
                    seed_dts.append(dt)
                    self.histograms["drift_time_seed"].Fill(dt_ns)
                    # 後の処理でパーセンタイルごとにFillするためのバッファ(inPixel用)は別途計算
                else:
                    neighbor_dts.append(dt)
                    self.histograms["drift_time_neighbor"].Fill(dt_ns)

            if all_dts:
                all_dts.sort()
                for c in self.target_carrier_counts:
                    if len(all_dts) >= c:
                        t_ns = all_dts[c - 1] * 1e-3
                        buffer[f"driftTime_to_{c}e"].append(t_ns)
                        buffer[f"inPixel_driftTime_to_{c}e"].append((in_x, in_y, t_ns))

            # 4. Percentile Calculation & Buffer Fill
            for p in self.percentiles:
                idx_func = lambda l, pct: min(int(len(l) * (pct/100.0)), len(l)-1)
                
                if all_dts:
                    all_dts.sort()
                    val = all_dts[idx_func(all_dts, p)] * 1e-3
                    buffer[f"electron_driftTime_{p}p"].append(val)
                    buffer[f"inPixel_electron_driftTime_{p}p"].append((in_x, in_y, val))
            
                if seed_dts:
                    seed_dts.sort()
                    val = seed_dts[idx_func(seed_dts, p)] * 1e-3
                    buffer[f"inPixel_driftTime_seed_{p}p"].append((in_x, in_y, val))
                
                if neighbor_dts:
                    neighbor_dts.sort()
                    val = neighbor_dts[idx_func(neighbor_dts, p)] * 1e-3
                    buffer[f"inPixel_driftTime_neighbor_{p}p"].append((in_x, in_y, val))

            # for p_obj in prop_objs:
            #     if p_obj.getType() != 255: continue # Electrons
            #     mcp_ptr = p_obj.getMCParticle()
            #     if not mcp_ptr or mcp_ptr.GetUniqueID() != target_mcp.particle_id: continue
                
            #     dt_ns = p_obj.getGlobalTime() * 1e-3
            #     mcp_id = mcp_ptr.GetUniqueID()
            #     dts.append(p_obj.getGlobalTime())

            #     if mcp_id in seed_mcp_ids:
            #         self.histograms["drift_time_seed"].Fill(dt_ns)
            #     else:
            #         self.histograms["drift_time_neighbor"].Fill(dt_ns)

                
            #     # ヒストグラムに直接Fill (bufferを経由しないほうがメモリに優しい)
            #     self.histograms["drift_time_spectrum"].Fill(dt_ns)
            #     if self.fill_3d:
            #         p_pos = p_obj.getLocalPosition()
            #         self.histograms["drift_time_map_xyz"].Fill(p_pos.X()*1000 - pix_center[0], 
            #                                                   p_pos.Y()*1000 - pix_center[1], 
            #                                                   p_pos.Z()*1000 - pix_center[2], dt_ns)

            # if dts:
            #     dts.sort()
            #     # t90 = dts[int(len(dts)*0.9)] * 1e-3
            #     # buffer["electron_driftTime_90p"].append(t90)
            #     # buffer["inPixel_electron_driftTime_90p"].append((in_x, in_y, t90))
            #     for p in self.percentiles:
            #         idx = min(int(len(dts) * (p / 100.0)), len(dts) - 1)
            #         tp_ns = dts[idx] * 1e-3
            #         buffer[f"electron_driftTime_{p}p"].append(tp_ns)
            #         buffer[f"inPixel_electron_driftTime_{p}p"].append((in_x, in_y, tp_ns))

            # --- 5. 強制クリーンアップ ---
            if (i_ev + 1) % BATCH_SIZE == 0:
                self._fill_histograms_from_buffer(buffer)
            
            # 各イベントの最後で大きなオブジェクトを削除
            del hits, mc_particles, primaries, clusters, dts
            
            if (i_ev + 1) % events_per_chunk == 0:
                gc.collect() 

        self._fill_histograms_from_buffer(buffer)
        gc.collect()

    # def run_analysis(self, events_per_chunk=5000):
    #     """
    #     Process events in smaller chunks to prevent memory buildup.
    #     """
    #     self.setup_histograms()
    #     branch = self.config.get("detector_name", "CE65")
    #     N_X, N_Y, MARGIN = 48, 24, 2
        
    #     # Memory management: fill histograms every BATCH_SIZE events
    #     # And trigger GC every events_per_chunk
    #     BATCH_SIZE = 100
    #     buffer = {k: [] for k in self.histograms}
        
    #     pitch = self.detector_model.pixel_size
    #     offset_vec = np.array([pitch[0]/2, pitch[1]/2, 0.0])

    #     print(f"Worker {self.worker_id}: Processing from {self.start_entry} to {self.end_entry}")

    #     for i_ev in tqdm(range(self.start_entry, self.end_entry), position=self.worker_id, leave=False):
    #         self.pixel_tree.GetEntry(i_ev)
    #         self.counters["Total Events"] += 1

    #         # (Analysis logic remains the same inside the loop...)
    #         # --- 1. MCParticles ---
    #         mcp_objs = getattr(self.mcp_tree, branch)
    #         mc_particles = {}
    #         primaries = []
    #         for obj in mcp_objs:
    #             pid = obj.GetUniqueID()
    #             ref = obj.getLocalReferencePoint()
    #             pos_um = np.array([ref.X(), ref.Y(), ref.Z()]) * 1000.0
    #             mcp = MCParticle(pid, obj.getParent().GetUniqueID() if obj.getParent() else 0, pos_um)
    #             mc_particles[pid] = mcp
    #             if mcp.parent == 0: primaries.append(mcp)

    #         if len(primaries) != 1:
    #             self.counters["Skipped: Multiple primary particles"] += 1
    #             continue
            
    #         target_mcp = primaries[0]
    #         particle_pos = target_mcp.local_reference_point + offset_vec
    #         ix_mc, iy_mc = self.detector_model.get_pixel_index(particle_pos)
            
    #         is_in_roi = (MARGIN <= ix_mc < N_X - MARGIN and MARGIN <= iy_mc < N_Y - MARGIN)
    #         in_pixel_pos_mc = particle_pos - self.detector_model.get_pixel_center(ix_mc, iy_mc)

    #         # --- 2. PixelHits & Clustering ---
    #         pixel_objs = getattr(self.pixel_tree, branch)
    #         hits = [PixelHit(o.getSignal(), o.getPixel().getIndex().X(), o.getPixel().getIndex().Y(),
    #                         [mc_particles[p.GetUniqueID()] for p in o.getMCParticles() if p.GetUniqueID() in mc_particles])
    #                 for o in pixel_objs]

    #         clusters = self.do_clustering(hits)
    #         best_cl = max(clusters, key=lambda c: c.seed_pixel_hit.signal) if clusters else None
            
    #         is_rec = False
    #         res_vec = None
    #         if best_cl:
    #             if target_mcp.particle_id in {p.particle_id for p in best_cl.get_mc_particles()}:
    #                 res_vec = particle_pos - best_cl.get_position(self.detector_model, self.one_bit_processing)
    #                 if np.linalg.norm(res_vec[:2]) <= 40: is_rec = True

    #         if is_in_roi:
    #             buffer["inPixel_efficiency"].append((in_pixel_pos_mc[0], in_pixel_pos_mc[1], 1.0 if is_rec else 0.0))

    #         if not is_rec:
    #             if not pixel_objs: self.counters["Skipped: No pixel hits"] += 1
    #             elif best_cl is None: pass
    #             elif res_vec is None: self.counters["Skipped: Cluster has no primary particle"] += 1
    #             else: self.counters["Skipped: Residual > 40 um"] += 1
    #             continue

    #         # --- 3. Accepted Cluster Detail Analysis ---
    #         seed = best_cl.seed_pixel_hit
    #         if not (MARGIN <= seed.pixel_index_x < N_X - MARGIN and MARGIN <= seed.pixel_index_y < N_Y - MARGIN):
    #             self.counters["Skipped: Hit at EDGE events"] += 1
    #             continue

    #         self.counters["Clusters Accepted"] += 1
    #         pix_center = self.detector_model.get_pixel_center(seed.pixel_index_x, seed.pixel_index_y)
    #         in_pix = particle_pos - pix_center
            
    #         q_sum = best_cl.charge
    #         q_seed = seed.signal
    #         q_nb_sum = q_sum - q_seed
    #         cl_size = best_cl.size
    #         res_r = np.linalg.norm(res_vec[:2])

    #         # Fill Buffer
    #         in_x, in_y = in_pix[0], in_pix[1]
    #         buffer["cluster_charge"].append(q_sum)
    #         buffer["cluster_size"].append(cl_size)
    #         buffer["seed_charge"].append(q_seed)
    #         buffer["residual_x"].append(res_vec[0])
    #         buffer["residual_y"].append(res_vec[1])
    #         buffer["residual_r"].append(res_r)
    #         buffer["inPixel_cluster_size"].append((in_x, in_y, cl_size))
    #         buffer["inPixel_seed_charge"].append((in_x, in_y, q_seed))
    #         buffer["inPixel_cluster_charge"].append((in_x, in_y, q_sum))
    #         buffer["inPixel_neighbor_charge_sum"].append((in_x, in_y, q_nb_sum))
    #         buffer["inPixel_residual_r"].append((in_x, in_y, res_r))

    #         if q_sum > 0:
    #             buffer["inPixel_seed_ratio"].append((in_x, in_y, q_seed / q_sum))
    #         buffer["inPixel_multi_hit_prob"].append((in_x, in_y, 1.0 if cl_size > 1 else 0.0))

    #         # Cluster Size Specific
    #         s_suffix = f"_size_{cl_size}" if cl_size <= self.max_cluster_size_hist else f"_size_{self.max_cluster_size_hist + 1}_plus"
    #         buffer[f"seed_charge{s_suffix}"].append(q_seed)
    #         buffer[f"residual_x{s_suffix}"].append(res_vec[0])
    #         buffer[f"residual_y{s_suffix}"].append(res_vec[1])
    #         buffer[f"residual_r{s_suffix}"].append(res_r)

    #         # --- 4. Drift Time ---
    #         prop_objs = getattr(self.propagated_tree, branch, [])
    #         dts = []
    #         for p_obj in prop_objs:
    #             if p_obj.getType() != 255: continue
    #             mcp_ptr = p_obj.getMCParticle()
    #             if not mcp_ptr or mcp_ptr.GetUniqueID() != target_mcp.particle_id: continue
                
    #             dt_ns = p_obj.getGlobalTime() * 1e-3
    #             dts.append(p_obj.getGlobalTime())
                
    #             p_pos = p_obj.getLocalPosition()
    #             create_um = np.array([p_pos.X(), p_pos.Y(), p_pos.Z()]) * 1000.0
    #             dist_um = np.linalg.norm(create_um - pix_center)
                
    #             self.histograms["drift_time_spectrum"].Fill(dt_ns)
    #             if self.fill_3d:
    #                 rel_pos = create_um - pix_center
    #                 self.histograms["drift_time_map_xyz"].Fill(rel_pos[0], rel_pos[1], rel_pos[2], dt_ns)

    #         if dts:
    #             dts.sort()
    #             t90 = dts[int(len(dts)*0.9)] * 1e-3
    #             buffer["electron_driftTime_90p"].append(t90)
    #             buffer["inPixel_electron_driftTime_90p"].append((in_x, in_y, t90))

    #         # --- Memory Release Logic ---
    #         if (i_ev + 1) % BATCH_SIZE == 0:
    #             self._fill_histograms_from_buffer(buffer)
            
    #         if (i_ev + 1) % events_per_chunk == 0:
    #             gc.collect() # Periodically clear memory within the worker process

    #     # Final fill for the remaining data
    #     self._fill_histograms_from_buffer(buffer)
    #     gc.collect()

    # def run_analysis(self):
    #     self.setup_histograms()
    #     branch = self.config.get("detector_name", "CE65")
    #     N_X, N_Y, MARGIN = 48, 24, 2
    #     BATCH_SIZE = 100
    #     buffer = {k: [] for k in self.histograms}
        
    #     pitch = self.detector_model.pixel_size
    #     offset_vec = np.array([pitch[0]/2, pitch[1]/2, 0.0])

    #     for i_ev in tqdm(range(self.start_entry, self.end_entry), position=self.worker_id, leave=False):
    #         self.pixel_tree.GetEntry(i_ev)
    #         self.counters["Total Events"] += 1

    #         # --- 1. MCParticles ---
    #         mcp_objs = getattr(self.mcp_tree, branch)
    #         mc_particles = {}
    #         primaries = []
    #         for obj in mcp_objs:
    #             pid = obj.GetUniqueID()
    #             ref = obj.getLocalReferencePoint()
    #             pos_um = np.array([ref.X(), ref.Y(), ref.Z()]) * 1000.0
    #             mcp = MCParticle(pid, obj.getParent().GetUniqueID() if obj.getParent() else 0, pos_um)
    #             mc_particles[pid] = mcp
    #             if mcp.parent == 0: primaries.append(mcp)

    #         if len(primaries) != 1:
    #             self.counters["Skipped: Multiple primary particles"] += 1
    #             continue
            
    #         target_mcp = primaries[0]
    #         particle_pos = target_mcp.local_reference_point + offset_vec
    #         ix_mc, iy_mc = self.detector_model.get_pixel_index(particle_pos)
            
    #         is_in_roi = (MARGIN <= ix_mc < N_X - MARGIN and MARGIN <= iy_mc < N_Y - MARGIN)
    #         in_pixel_pos_mc = particle_pos - self.detector_model.get_pixel_center(ix_mc, iy_mc)

    #         # --- 2. PixelHits & Clustering ---
    #         pixel_objs = getattr(self.pixel_tree, branch)
    #         hits = [PixelHit(o.getSignal(), o.getPixel().getIndex().X(), o.getPixel().getIndex().Y(),
    #                         [mc_particles[p.GetUniqueID()] for p in o.getMCParticles() if p.GetUniqueID() in mc_particles])
    #                 for o in pixel_objs]

    #         clusters = self.do_clustering(hits)
    #         best_cl = max(clusters, key=lambda c: c.seed_pixel_hit.signal) if clusters else None
            
    #         is_rec = False
    #         res_vec = None
    #         if best_cl:
    #             if target_mcp.particle_id in {p.particle_id for p in best_cl.get_mc_particles()}:
    #                 res_vec = particle_pos - best_cl.get_position(self.detector_model, self.one_bit_processing)
    #                 if np.linalg.norm(res_vec[:2]) <= 40: is_rec = True

    #         if is_in_roi:
    #             buffer["inPixel_efficiency"].append((in_pixel_pos_mc[0], in_pixel_pos_mc[1], 1.0 if is_rec else 0.0))

    #         if not is_rec:
    #             if not pixel_objs: self.counters["Skipped: No pixel hits"] += 1
    #             elif best_cl is None: pass
    #             elif res_vec is None: self.counters["Skipped: Cluster has no primary particle"] += 1
    #             else: self.counters["Skipped: Residual > 40 um"] += 1
    #             continue

    #         # --- 3. Accepted Cluster Detail Analysis ---
    #         seed = best_cl.seed_pixel_hit
    #         if not (MARGIN <= seed.pixel_index_x < N_X - MARGIN and MARGIN <= seed.pixel_index_y < N_Y - MARGIN):
    #             self.counters["Skipped: Hit at EDGE events"] += 1
    #             continue

    #         self.counters["Clusters Accepted"] += 1
    #         pix_center = self.detector_model.get_pixel_center(seed.pixel_index_x, seed.pixel_index_y)
    #         in_pix = particle_pos - pix_center
            
    #         q_sum = best_cl.charge
    #         q_seed = seed.signal
    #         q_nb_sum = q_sum - q_seed
    #         cl_size = best_cl.size
    #         res_r = np.linalg.norm(res_vec[:2])

    #         # Standard Histograms
    #         buffer["cluster_charge"].append(q_sum)
    #         buffer["cluster_size"].append(cl_size)
    #         buffer["seed_charge"].append(q_seed)
    #         buffer["residual_x"].append(res_vec[0])
    #         buffer["residual_y"].append(res_vec[1])
    #         buffer["residual_r"].append(res_r)
            
    #         # In-Pixel Profile fills
    #         in_x, in_y = in_pix[0], in_pix[1]
    #         buffer["inPixel_cluster_size"].append((in_x, in_y, cl_size))
    #         buffer["inPixel_seed_charge"].append((in_x, in_y, q_seed))
    #         buffer["inPixel_cluster_charge"].append((in_x, in_y, q_sum))
    #         buffer["inPixel_neighbor_charge_sum"].append((in_x, in_y, q_nb_sum))
    #         buffer["inPixel_residual_r"].append((in_x, in_y, res_r))
    #         buffer["inPixel_residual_x"].append((in_x, in_y, abs(res_vec[0])))
    #         buffer["inPixel_residual_y"].append((in_x, in_y, abs(res_vec[1])))

    #         if q_sum > 0:
    #             buffer["inPixel_seed_ratio"].append((in_x, in_y, q_seed / q_sum))
            
    #         buffer["inPixel_multi_hit_prob"].append((in_x, in_y, 1.0 if cl_size > 1 else 0.0))

    #         # Cluster Size Specific Histograms
    #         if cl_size <= self.max_cluster_size_hist:
    #             s_suffix = f"_size_{cl_size}"
    #         else:
    #             s_suffix = f"_size_{self.max_cluster_size_hist + 1}_plus"
            
    #         buffer[f"seed_charge{s_suffix}"].append(q_seed)
    #         buffer[f"residual_x{s_suffix}"].append(res_vec[0])
    #         buffer[f"residual_y{s_suffix}"].append(res_vec[1])
    #         buffer[f"residual_r{s_suffix}"].append(res_r)

    #         # --- 4. Propagated Charge (Drift Time) ---
    #         prop_objs = getattr(self.propagated_tree, branch, [])
    #         dts = []
    #         for p_obj in prop_objs:
    #             if p_obj.getType() != 255: continue # Electrons only
                
    #             # Check link to primary particle
    #             mcp_ptr = p_obj.getMCParticle()
    #             if not mcp_ptr or mcp_ptr.GetUniqueID() != target_mcp.particle_id:
    #                 continue
                
    #             dt_ns = p_obj.getGlobalTime() * 1e-3
    #             dts.append(p_obj.getGlobalTime()) # ps
                
    #             p_pos = p_obj.getLocalPosition()
    #             create_um = np.array([p_pos.X(), p_pos.Y(), p_pos.Z()]) * 1000.0
    #             dist_um = np.linalg.norm(create_um - pix_center)
                
    #             # Direct fill for drift histograms to save memory
    #             if self.histograms.get("drift_time_vs_distance"): 
    #                 self.histograms["drift_time_vs_distance"].Fill(dist_um, dt_ns)
    #             if self.histograms.get("drift_time_spectrum"):
    #                 self.histograms["drift_time_spectrum"].Fill(dt_ns)
    #             if self.histograms.get("drift_time_vs_depth"):
    #                 self.histograms["drift_time_vs_depth"].Fill(create_um[2], dt_ns)

    #             if self.fill_3d and self.histograms.get("drift_time_map_xyz"):
    #                 rel_pos = create_um - pix_center
    #                 self.histograms["drift_time_map_xyz"].Fill(rel_pos[0], rel_pos[1], rel_pos[2], dt_ns)

    #         if dts:
    #             dts.sort()
    #             t90 = dts[int(len(dts)*0.9)] * 1e-3 # 90% drift time in ns
    #             buffer["electron_driftTime_90p"].append(t90)
    #             buffer["inPixel_electron_driftTime_90p"].append((in_x, in_y, t90))
    #             buffer["cluster_charge_vs_drift_time"].append((q_sum, t90))
    #             buffer["cluster_size_vs_drift_time"].append((cl_size, t90))

    #         if (i_ev + 1) % BATCH_SIZE == 0:
    #             self._fill_histograms_from_buffer(buffer)
    #             gc.collect()

    #     self._fill_histograms_from_buffer(buffer)

    def finalize(self):
        # for i, name in enumerate(self.counter_names):
        #     count = self.counters.get(name, 0)
        #     self.histograms["counters"].SetBinContent(i + 1, count)
        
        # output_filename = self.config.get("output_file_name", "analysis_py.root")
        # output_file = ROOT.TFile(output_filename, "RECREATE")
        # for name, histo in self.histograms.items():
        #     histo.Write()
        # output_file.Close()
        # print(f"Histograms written to {output_filename}")

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
            dir_name = f"clsize_{i}"
            dir_size_specific[f"_size_{i}"] = dir_per_size.mkdir(dir_name)

        suffix_plus = f"{self.max_cluster_size_hist + 1}_plus"
        dir_name_plus = f"clsize_{self.max_cluster_size_hist + 1}_plus"
        dir_size_specific[f"_size_{suffix_plus}"] = dir_per_size.mkdir(dir_name_plus)

        for name, histo in self.histograms.items():
            
            written_to_subdir = False
            
            # --- 分類ルール1: クラスターサイズ別 ---
            for suffix, directory in dir_size_specific.items():
                if name.endswith(suffix):
                    directory.cd()      # サブディレクトリに移動
                    histo.Write()
                    written_to_subdir = True
                    break # このヒストグラムの処理は完了
            
            if written_to_subdir:
                continue # 次のヒストグラムへ

            # --- 分類ルール2: inPixel ---
            if name.startswith("inPixel_"):
                dir_inpixel.cd()
                histo.Write()
                continue

            # --- 分類ルール3: DriftTime ---
            if (name.startswith("drift_time_") or 
                name.startswith("electron_driftTime_") or 
                name.endswith("_vs_drift_time") or
                name.endswith("_vs_depth")):
                dir_drifttime.cd()
                histo.Write()
                continue

            # --- 分類ルール4: 2D相関 ---
            # if "_vs_" in name:
            #     dir_2d_correlations.cd()
            #     histo.Write()
            #     continue

            # --- 分類ルール5: その他 (ルートディレクトリ) ---
            output_file.cd() # メインディレクトリに戻る
            histo.Write()

        # --- 4. ファイルを閉じる ---
        output_file.Close()
        print(f"Histograms written to {output_filename}")

def run_worker(args):
    start_index, end_index, worker_id, config = args

    # ワーカーごとに固有の出力ファイル名とイベント範囲を設定
    config["output_file_name"] = f"analysis_py_part_{worker_id}.root"
    config["start_entry"] = start_index
    config["end_entry"] = end_index
    config["worker_id"] = worker_id # デバッグ用

    try:
        analyzer = AnalysisPixelModule(config)
        analyzer.run_analysis(events_per_chunk=2000)
        analyzer.finalize()
        return config["output_file_name"] # 成功したらファイル名を返す
    except Exception as e:
        # エラーが発生した場合、トレースバックを表示
        print(f"--- Worker {worker_id} FAILED (Events {start_index}-{end_index}) ---")
        print(f"An error occurred in worker {worker_id}: {e}")
        traceback.print_exc()
        print(f"---------------------------------")
        return None # 失敗したらNoneを返す

if __name__ == '__main__':
    start_time = time.time()
    
    parser = argparse.ArgumentParser(
        description="Allpix Squared Analysis script (Multiprocessing Enabled).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 引数のマージ ---
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to the input ROOT file from Allpix Squared Simulation Data.")
    parser.add_argument("-o", "--output", default="analysis_py.root", type=str, help="Path for the FINAL merged output ROOT file.")
    parser.add_argument("-st", "--seed_threshold", type=int, default=0, help="Seed_threshold for analysis (electron num).")
    parser.add_argument("-nt", "--neighbor_threshold", type=int, default=0, help="Neighbor_threshold for analysis (electron num).")
    parser.add_argument("-p", "--pixel_pitch", type=float, default=22.5, help="pixel pitch about DUT")
    parser.add_argument("-v", "--voltage", type=int, default=10, help="Chip voltage about DUT (e.g., 10, 7, 4)")
    parser.add_argument("-b", "--beam_type", type=str, default="e3GeV", help="Beam information e.g., e3GeV")
    parser.add_argument("-m", "--model", type=str, default="masetti", help="Model name for Electron calculation")
    parser.add_argument("-n", "--name", type=str, default="CE65", help="DUT name (TTree branch name)")
    
    # ベースコードの引数
    parser.add_argument("--one_bit", action="store_true", help="Enable 1-bit processing (binary readout)")
    
    # 追加機能コードの引数
    parser.add_argument("--fill3D", action="store_true", help="Enable filling of the 3D drift time map (can be slow).")
    parser.add_argument("-j", "--cores", type=int, default=8, help="Number of CPU cores to use for multiprocessing.")

    args = parser.parse_args()    

    # --- config辞書のマージ ---
    config = {
        "file_name": args.input_file,
        "output_file_name": args.output, # これは最終ファイル名として使う
        "detector_name": args.name,
        "pixel_pitch_x": args.pixel_pitch,
        "pixel_pitch_y": args.pixel_pitch,
        "seed_threshold": args.seed_threshold,
        "neighbor_threshold": args.neighbor_threshold,
        "voltage": args.voltage,
        "beam_type": args.beam_type,
        "model_name": args.model,
        "granularity_x": 50,
        "granularity_y": 50, 
        "max_cluster_charge": 20000, # 60keに変更
        "one_bit": args.one_bit,     # ベースコードから追加
        "fill_3d": args.fill3D,    # 追加機能コードから追加
        "percentiles": [50, 60, 70, 80, 90],
        "target_carrier_counts": [48, 96, 144, 150, 160, 170, 140, 180, 190,200, 240, 312],
    }

    num_cores = args.cores

    # --- タスク分割 ---
    print(f"Opening file {config['file_name']} to get total entries...")
    try:
        temp_file = ROOT.TFile.Open(config["file_name"])
        if not temp_file or temp_file.IsZombie():
            raise RuntimeError(f"File not found or is zombie: {config['file_name']}")
        
        # 必要なTTreeが存在するかチェック
        if not temp_file.Get("PixelHit") or not temp_file.Get("MCParticle") or not temp_file.Get("PropagatedCharge"):
             print(f"Error: Input file '{config['file_name']}' is missing required TTrees.")
             print(f"Required: 'PixelHit', 'MCParticle', AND 'PropagatedCharge'")
             temp_file.Close()
             exit(1)
             
        temp_tree = temp_file.Get("PixelHit")
        total_entries = temp_tree.GetEntries()
        temp_file.Close()
    except Exception as e:
        print(f"Error opening file to get entries: {e}")
        exit(1)

    if total_entries == 0:
        print("Error: Input file has 0 entries. Exiting.")
        exit(0)

    print(f"Total entries: {total_entries}. Splitting tasks for {num_cores} cores.")

    chunk_size = math.ceil(total_entries / num_cores)
    tasks = [] # (start_index, end_index, worker_id, config_copy)

    for i in range(num_cores):
        start_index = i * chunk_size
        end_index = min((i + 1) * chunk_size, total_entries)
        if start_index >= total_entries:
            continue # イベントがもう残っていない
        
        # 各ワーカーに渡す引数を準備 (configはコピーして渡す)
        tasks.append((start_index, end_index, i, config.copy()))

    if not tasks:
        print("No tasks to run. Exiting.")
        exit(0)

    # --- プールを開始してタスクを実行 ---
    print(f"Starting multiprocessing pool with {len(tasks)} workers...")
    
    partial_files = [] # 成功したワーカーが出力したファイル名のリスト
    
    # .Pool()のコンテキストマネージャを使用
    with multiprocessing.Pool(processes=num_cores) as pool:
        # pool.map が全ワーカーの終了を待機し、結果 (ファイル名 or None) のリストを返す
        results = pool.map(run_worker, tasks)
        
        for res in results:
            if res: # Noneでない（＝成功した）場合のみリストに追加
                partial_files.append(res)

    if not partial_files:
        print("All workers failed. No partial files were produced. Exiting.")
        exit(1)
        
    print(f"\nAll workers finished. Found {len(partial_files)} partial files.")

    # --- haddによるマージ ---
    final_output_file = config["output_file_name"]
    hadd_command = f"hadd -f {final_output_file} {' '.join(partial_files)}"
    
    print("Merging partial files with hadd...")
    print(f"> {hadd_command}")
    
    try:
        os.system(hadd_command)
        print(f"Successfully merged into {final_output_file}")
    except Exception as e:
        print(f"hadd command FAILED: {e}")
        print("Partial files are left for inspection.")
        exit(1)

    # --- クリーンアップ: 部分ファイルを削除 ---
    print("Cleaning up partial files...")
    for f in partial_files:
        try:
            os.remove(f)
        except OSError as e:
            print(f"Warning: could not remove partial file {f}: {e}")

    # --- 実行時間レポート ---
    end_time = time.time()
    elapsed_time = end_time - start_time
    start_str = datetime.fromtimestamp(start_time).strftime("%Y/%m/%d %H:%M:%S")
    end_str = datetime.fromtimestamp(end_time).strftime("%Y/%m/%d %H:%M:%S")
    print("------------------------------------------")
    print(f"Start time : {start_str}")
    print(f"End time   : {end_str}")
    print(f"Total time : {elapsed_time:.2f} seconds (including merge)")

    minutes = int(elapsed_time // 60)
    seconds = elapsed_time % 60
    if minutes > 0:
        print(f"Which is:    {minutes} minute(s) and {seconds:.2f} second(s)")
    print("------------------------------------------")