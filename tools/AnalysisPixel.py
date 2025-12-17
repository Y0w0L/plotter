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
        
        # particle_idを使って重複を削除し、順番を保持する
        unique_particles = {}
        for particle in all_particles:
            if particle.particle_id not in unique_particles:
                unique_particles[particle.particle_id] = particle
        
        return list(unique_particles.values())


class DetectorModel:
    """検出器のジオメトリ情報を保持するヘルパークラス"""
    def __init__(self, pixel_size_x_um, pixel_size_y_um):
        self.pixel_size = np.array([pixel_size_x_um, pixel_size_y_um])

    def get_pixel_center(self, ix, iy):
        return np.array([(ix + 0.5) * self.pixel_size[0], (iy + 0.5) * self.pixel_size[1], 0.0])

    def are_neighbors(self, index1_x, index1_y, index2_x, index2_y, distance=1):
        return abs(index1_x - index2_x) <= distance and abs(index1_y - index2_y) <= distance

    def get_pixel_index(self, position):
        #ix = int(position[0] / self.pixel_size[0])
        #iy = int(position[1] / self.pixel_size[1])
        #return ix, iy
        ix = math.floor(position[0] / self.pixel_size[0])
        iy = math.floor(position[1] / self.pixel_size[1])
        return ix, iy

class AnalysisPixelModule:
    """C++のAnalysisPixelModuleのロジックを再現するPythonクラス"""
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

        # PyROOTでファイルとTTreeを開く
        self.input_file = ROOT.TFile.Open(config["file_name"])
        self.pixel_tree = self.input_file.Get("PixelHit")
        self.mcp_tree = self.input_file.Get("MCParticle")
        self.n_entries = self.pixel_tree.GetEntries()
        #self.n_entries = 100  
        self.propagated_tree = self.input_file.Get("PropagatedCharge")

        if not self.pixel_tree or not self.mcp_tree or not self.propagated_tree:
            raise RuntimeError(f"File '{config['file_name']}' is missing required TTrees (PixelHit, MCParticle, or PropagatedCharge).")

        self.n_entries_total = self.pixel_tree.GetEntries()
        self.start_entry = config.get("start_entry", 0)
        # configでend_entryが指定されなければ、全イベントを処理
        self.end_entry = config.get("end_entry", self.n_entries_total)
        # このワーカーが処理するイベント数
        self.n_entries = self.end_entry - self.start_entry

        # 2つのTTreeが同じイベントを指すように同期させる
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

        #self.four_neighbors = [(0, 1), (0, -1), (1, 0), (-1, 0)]
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
        self.histograms["inPixel_efficiency"] = ROOT.TProfile2D("inPixel_efficiency", ";x/pitch [um];y/pitch [um];efficiency", inpixel_bins_x, -pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
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

        self.histograms["prof_clusterSize_vs_clusterCharge"] = ROOT.TProfile(
            "prof_clusterSize_vs_clusterCharge",
            "Average Cluster Size vs Cluster Charge;Cluster Charge [e];Average Cluster Size [pixels]",
            100, 0, max_cluster_charge
        )

        self.histograms["inPixel_seed_ratio"] = ROOT.TProfile2D(
            "inPixel_seed_ratio", 
            ";x/pitch [um];y/pitch [um];Seed Charge Ratio (Q_seed / Q_clus)", 
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )
        
        self.histograms["inPixel_multi_hit_prob"] = ROOT.TProfile2D(
            "inPixel_multi_hit_prob", 
            ";x/pitch [um];y/pitch [um];Multi-Pixel Hit Probability", 
            inpixel_bins_x, -pitch_x / 2, pitch_x / 2, 
            inpixel_bins_y, -pitch_y / 2, pitch_y / 2
        )

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

        self.histograms["inPixel_electron_driftTime_90p"] = ROOT.TProfile2D("inPixel_electron_driftTime_90p", ";x/pitch [um];y/pitch [um];90% electron drift time [ns]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["electron_driftTime_90p"] = ROOT.TH1D("electron_driftTime_90p", ";90% electron drift time [ns];counts", 20000, 0, 2)

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

        self.histograms["drift_time_spectrum"] = ROOT.TH1D("drift_time_spectrum", "All Electron Drift Time;Drift Time [ns];Counts", 2000, 0, 0.15)
        self.histograms["drift_time_vs_depth"] = ROOT.TH2D("drift_time_vs_depth", "Drift Time vs Creation Depth;Creation Z [um];drift time [ns]", 100, -30, 30, 2000, 0, 0.15)

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
        # (この関数は変更ありません)
        clusters = []
        used_pixels = set()

        hit_map = {(p.pixel_index_x, p.pixel_index_y): p for p in pixel_hits} # add for acceleration

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
                # for neighbor_candidate in sorted_pixel_hits:
                #     neighbor_id = (neighbor_candidate.pixel_index_x, neighbor_candidate.pixel_index_y)
                #     if neighbor_id in used_pixels or neighbor_candidate.signal < self.neighbor_threshold:
                #         continue
                #     if self.detector_model.are_neighbors(current_pixel.pixel_index_x, current_pixel.pixel_index_y, neighbor_candidate.pixel_index_x, neighbor_candidate.pixel_index_y):
                #         cluster.add_pixel_hit(neighbor_candidate)
                #         used_pixels.add(neighbor_id)
                #         to_check.append(neighbor_candidate)

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
        # (この関数は変更ありません)
        return [p for p in mc_particles.values() if p.parent is None]

    def _fill_histograms_from_buffer(self, buffer: Dict[str, list]):
        """バッファに溜まったデータからヒストグラムを更新する"""
        # for key in buffer:
        #     buffer[key].clear()
        for key, histo in self.histograms.items():
            if not buffer.get(key): continue
            
            # ヒストグラムの型に応じてFillメソッドを呼び分ける
            try:
                if isinstance(histo, ROOT.TProfile3D):
                    for x, y, z, val in buffer[key]:
                        histo.Fill(x, y, z, val)
                elif isinstance(histo, ROOT.TProfile2D):
                    for x, y, val in buffer[key]:
                        histo.Fill(x, y, val)
                elif isinstance(histo, ROOT.TH2D):
                    for x_val, y_val in buffer[key]:
                        histo.Fill(x_val, y_val)
                elif isinstance(histo, ROOT.TH1D):
                    for val in buffer[key]:
                        histo.Fill(val)
            except Exception as e:
                print(f"Error filling histogram {key}: {e}")
                # エラーが発生しても処理を続ける
        
        # バッファをクリア
        for key in buffer:
            buffer[key].clear()

    def run_analysis(self):
        """PyROOTを使ってイベントループを実行し、解析を行う"""
        self.setup_histograms()
        
        branch_name = "CE65"

        N_PIXELS_X = 48
        N_PIXELS_Y = 24
        MARGIN = 2

        self.counters = {name: 0 for name in self.counter_names}
        self.counters["Total Events"] = self.n_entries

        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        offset_x = pitch_x / 2
        offset_y = pitch_y / 2
        
        BATCH_SIZE = 1000
        buffer = {hist_name: [] for hist_name in self.histograms}

        # 高速化のために変数をキャッシュ
        h_drift_time_vs_dist = self.histograms.get("drift_time_vs_distance")
        h_drift_time_map_xyz = self.histograms.get("drift_time_map_xyz") if self.fill_3d else None
        h_drift_time_seed = self.histograms.get("drift_time_seed")
        h_drift_time_neighbor = self.histograms.get("drift_time_neighbor")
        h_drift_time_spectrum = self.histograms.get("drift_time_spectrum")
        h_drift_time_vs_depth = self.histograms.get("drift_time_vs_depth")

        tqdm_iterator = tqdm(
            range(self.start_entry, self.end_entry),
            desc=f"Worker {self.worker_id} (Events {self.start_entry}-{self.end_entry})", 
            position=self.worker_id,
            leave=False
        )

        for i_global in tqdm_iterator:
            self.pixel_tree.GetEntry(i_global)

            # --- 1. MCParticle データデコード (位置を最初に特定するため移動) ---
            mcp_objects = getattr(self.mcp_tree, branch_name)
            mc_particles: Dict[int, MCParticle] = {}
            
            for mcp_obj in mcp_objects:
                pid = mcp_obj.GetUniqueID()
                parent_ptr = mcp_obj.getParent()
                parent_id = parent_ptr.GetUniqueID() if parent_ptr else 0
                ref_point = mcp_obj.getLocalReferencePoint()
                mc_particles[pid] = MCParticle(
                    particle_id=pid,
                    parent=parent_id,
                    local_reference_point=np.array([ref_point.X()*1000, ref_point.Y()*1000, ref_point.Z()*1000])
                )

            for mcp in mc_particles.values():
                parent_id = mcp.parent
                if parent_id != 0 and parent_id in mc_particles:
                    mcp.parent = mc_particles[parent_id]
                else:
                    mcp.parent = None

            # --- 2. Primary Particle チェックとROI判定 ---
            primary_particles = self.get_primary_particles(mc_particles)
            
            if len(primary_particles) != 1:
                self.counters["Skipped: Multiple primary particles"] += 1
                continue
            
            # Primary粒子の位置を特定
            target_particle = primary_particles[0]
            particle_pos = target_particle.local_reference_point + np.array([offset_x, offset_y, 0])
            
            # MC粒子がROI内にあるか確認 (Efficiencyの分母)
            ix_mc, iy_mc = self.detector_model.get_pixel_index(particle_pos)
            
            is_in_roi = (MARGIN <= ix_mc < N_PIXELS_X - MARGIN and MARGIN <= iy_mc < N_PIXELS_Y - MARGIN)
            
            pixel_center_mc = self.detector_model.get_pixel_center(ix_mc, iy_mc)
            in_pixel_pos_mc = particle_pos - pixel_center_mc

            # --- 3. PixelHit 処理とクラスタリング ---
            pixel_objects = getattr(self.pixel_tree, branch_name)
            
            # Efficiency計算用のフラグ
            is_reconstructed_successfully = False
            best_cluster = None
            residual_vec = None
            residual_r = None

            # ヒットがある場合のみクラスタリング処理を行う
            if len(pixel_objects) > 0:
                pixel_hits = []
                for hit_obj in pixel_objects:
                    hit_mc_particles = []
                    for linked_mcp_ptr in hit_obj.getMCParticles():
                        linked_pid = linked_mcp_ptr.GetUniqueID()
                        if linked_pid in mc_particles:
                            hit_mc_particles.append(mc_particles[linked_pid])
                    
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
                    
                    # クラスタがPrimary粒子と紐付いているか確認
                    cluster_particles_ids = {p.particle_id for p in best_cluster.get_mc_particles()}
                    if target_particle.particle_id in cluster_particles_ids:
                        
                        # 残差計算
                        cluster_pos = best_cluster.get_position(self.detector_model, one_bit=self.one_bit_processing)
                        residual_vec = particle_pos - cluster_pos
                        residual_r = np.linalg.norm(residual_vec[:2])

                        cut_um = 40
                        if residual_r <= cut_um:
                            is_reconstructed_successfully = True

            # --- 4. Efficiency Plot の Fill ---
            # ROI内に入射したイベントであれば、検出できた(1)か否(0)かを記録
            if is_in_roi:
                eff_val = 1.0 if is_reconstructed_successfully else 0.0
                buffer["inPixel_efficiency"].append((in_pixel_pos_mc[0], in_pixel_pos_mc[1], eff_val))

            # --- 5. 検出失敗ならここで終了 (Efficiency 0 のケース) ---
            if not is_reconstructed_successfully:
                if len(pixel_objects) == 0:
                    self.counters["Skipped: No pixel hits"] += 1
                elif best_cluster is None:
                    # ヒットはあるがクラスタにならなかった（閾値など）
                    pass 
                elif residual_vec is None:
                    # クラスタはあるがPrimaryと紐付かなかった
                    self.counters["Skipped: Cluster has no primary particle"] += 1
                else:
                    # 残差カット落ち
                    self.counters["Skipped: Residual > 40 um"] += 1
                continue

            # --- 6. 以下、検出成功 (Efficiency 1) イベントの詳細解析 ---
            # ここから先は `best_cluster` が有効で、かつResidualカットも通過していることが保証される
            
            # Seed PixelがEDGEにある場合は除く（従来の解析ロジックを踏襲）
            # 注意: EfficiencyMapはMC位置基準でROIを切ったが、詳細プロットはSeed位置基準で切るのが通例
            seed_x = best_cluster.seed_pixel_hit.pixel_index_x
            seed_y = best_cluster.seed_pixel_hit.pixel_index_y

            if not (MARGIN <= seed_x < N_PIXELS_X - MARGIN and MARGIN <= seed_y < N_PIXELS_Y - MARGIN):
                self.counters["Skipped: Hit at EDGE events"] += 1
                continue

            self.counters["Clusters Accepted"] += 1
            clus = best_cluster # 変数名を合わせる

            # Propagated Charge の処理 (Drift Time用)
            propagated_objects = getattr(self.propagated_tree, branch_name, [])
            charge_map: Dict[int, List[PropagatedCharge]] = {}
            for prop_obj in propagated_objects:
                mc_particle_ptr = prop_obj.getMCParticle()
                if not mc_particle_ptr: continue
                parent_mcp = mc_particles.get(mc_particle_ptr.GetUniqueID())
                if not parent_mcp: continue
                
                carrier_type_enum = prop_obj.getType()
                if carrier_type_enum == 255: charge_type = 'electron'
                elif carrier_type_enum == 1: charge_type = 'hole'
                else: continue

                pos = prop_obj.getLocalPosition()
                creation_pos_um = np.array([pos.X() * 1000, pos.Y() * 1000, pos.Z() * 1000])
                
                prop_charge = PropagatedCharge(
                    charge_type=charge_type,
                    global_time=prop_obj.getGlobalTime(),
                    mc_particle=parent_mcp,
                    local_creation_pos=creation_pos_um
                )
                charge_map.setdefault(parent_mcp.particle_id, []).append(prop_charge)


            # In-Pixel Positionの再計算 (Seed基準: プロットの一貫性のため)
            ix, iy = self.detector_model.get_pixel_index(particle_pos)
            pixel_center = self.detector_model.get_pixel_center(ix, iy)
            in_pixel_pos = particle_pos - pixel_center
            
            non_seed_hits = [hit for hit in clus.pixel_hits if hit is not clus.seed_pixel_hit]
            sum_non_seed_charge = sum(hit.signal for hit in non_seed_hits)

            # ヒストグラムへのFill
            buffer["cluster_charge"].append(clus.charge)
            buffer["cluster_size"].append(clus.size)
            buffer["seed_charge"].append(clus.seed_pixel_hit.signal)
            buffer["residual_x"].append(residual_vec[0])
            buffer["residual_y"].append(residual_vec[1])
            buffer["residual_r"].append(residual_r)
            buffer["inPixel_cluster_size"].append((in_pixel_pos[0], in_pixel_pos[1], clus.size))
            buffer["inPixel_residual_r"].append((in_pixel_pos[0], in_pixel_pos[1], residual_r))
            buffer["inPixel_residual_x"].append((in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[0])))
            buffer["inPixel_residual_y"].append((in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[1])))
            buffer["inPixel_residual_xy2"].append((in_pixel_pos[0], in_pixel_pos[1], (abs(residual_vec[0]) + abs(residual_vec[1])) / 2))
            buffer["inPixel_seed_charge"].append((in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal))
            buffer["inPixel_cluster_charge"].append((in_pixel_pos[0], in_pixel_pos[1], clus.charge))
            buffer["inPixel_neighbor_charge_sum"].append((in_pixel_pos[0], in_pixel_pos[1], sum_non_seed_charge))
            buffer["clusterSize_vs_clusterCharge"].append((clus.charge, clus.size))
            buffer["prof_clusterSize_vs_clusterCharge"].append((clus.charge, clus.size))
           
            # 1. Seed Charge Ratio (0.0 ~ 1.0)
            if clus.charge > 0:
                ratio = clus.seed_pixel_hit.signal / clus.charge
                buffer["inPixel_seed_ratio"].append((in_pixel_pos[0], in_pixel_pos[1], ratio))

            # 2. Multi-Pixel Hit Probability (0 or 1)
            # サイズが1より大きければ 1.0 (True), そうでなければ 0.0 (False)
            is_multi = 1.0 if clus.size > 1 else 0.0
            buffer["inPixel_multi_hit_prob"].append((in_pixel_pos[0], in_pixel_pos[1], is_multi))


            if sum_non_seed_charge > 0:
                buffer["cluster_neighbor_charge_sum"].append(sum_non_seed_charge)
                buffer["neighborChargeSum_vs_clusterSize"].append((sum_non_seed_charge, clus.size))
                buffer["seedCharge_vs_neighborChargeSum"].append((clus.seed_pixel_hit.signal, sum_non_seed_charge))

            buffer["seedCharge_vs_clusterSize"].append((clus.seed_pixel_hit.signal, clus.size))
            buffer["clusterCharge_vs_clusterSize"].append((clus.charge, clus.size))

            for hit in non_seed_hits:
                buffer["cluster_neighbor_charge"].append(hit.signal)

            size = clus.size
            seed_charge_ke = clus.seed_pixel_hit.signal

            if size <= self.max_cluster_size_hist:
                buffer[f"seed_charge_size_{size}"].append(seed_charge_ke)
                buffer[f"residual_x_size_{size}"].append(residual_vec[0])
                buffer[f"residual_y_size_{size}"].append(residual_vec[1])
                buffer[f"residual_r_size_{size}"].append(residual_r)
            else:
                suffix_plus = f"{self.max_cluster_size_hist + 1}_plus"
                buffer[f"seed_charge_size_{suffix_plus}"].append(seed_charge_ke)
                buffer[f"residual_x_size_{suffix_plus}"].append(residual_vec[0])
                buffer[f"residual_y_size_{suffix_plus}"].append(residual_vec[1])
                buffer[f"residual_r_size_{suffix_plus}"].append(residual_r)
            
            # --- Drift Time Analysis ---
            electron_drift_times_ps = []
            
            for hit in clus.pixel_hits:
                pixel_center_pos = self.detector_model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
                
                for mcp in hit.mc_particles:
                    propagated_charges = charge_map.get(mcp.particle_id, [])
                    
                    for pc in propagated_charges:
                        if pc.charge_type == 'electron':
                            drift_time_ns = pc.global_time * 1e-3
                            electron_drift_times_ps.append(pc.global_time)
                            
                            creation_pos = pc.local_creation_pos
                            distance_um = np.linalg.norm(creation_pos - pixel_center_pos)
                            
                            if h_drift_time_vs_dist: h_drift_time_vs_dist.Fill(distance_um, drift_time_ns)
                            if h_drift_time_spectrum: h_drift_time_spectrum.Fill(drift_time_ns)
                            if h_drift_time_vs_depth: h_drift_time_vs_depth.Fill(creation_pos[2], drift_time_ns)

                            if h_drift_time_map_xyz:
                                relative_pos = creation_pos - pixel_center_pos
                                h_drift_time_map_xyz.Fill(relative_pos[0], relative_pos[1], relative_pos[2], drift_time_ns)
                            
                            if hit is clus.seed_pixel_hit:
                                if h_drift_time_seed: h_drift_time_seed.Fill(drift_time_ns)
                            else:
                                if h_drift_time_neighbor: h_drift_time_neighbor.Fill(drift_time_ns)

            if electron_drift_times_ps:
                electron_drift_times_ps.sort()
                idx_90 = int(len(electron_drift_times_ps) * 0.9)
                time_90p_ns = electron_drift_times_ps[idx_90] * 1e-3

                buffer["electron_driftTime_90p"].append(time_90p_ns)
                buffer["inPixel_electron_driftTime_90p"].append((in_pixel_pos[0], in_pixel_pos[1], time_90p_ns))
                buffer["cluster_charge_vs_drift_time"].append((clus.charge, time_90p_ns))
                buffer["cluster_size_vs_drift_time"].append((clus.size, time_90p_ns))

            if (i_global + 1) % BATCH_SIZE == 0:
                self._fill_histograms_from_buffer(buffer)
        
        if any(buffer.values()):
            self._fill_histograms_from_buffer(buffer)

        print("Analysis finished.")
        print("----------------- Report -----------------")
        for key, value in self.counters.items():
            print(f"{key:<40}:{value}")
        print("------------------------------------------")
        

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
        analyzer.run_analysis()
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
    parser.add_argument("-j", "--cores", type=int, default=6, help="Number of CPU cores to use for multiprocessing.")

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