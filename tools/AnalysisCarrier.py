import ROOT
import numpy as np
import math
from dataclasses import dataclass, field
from typing import List, Dict, Any
from ROOT import gSystem
from tqdm import tqdm
import argparse
import time
from datetime import datetime

gSystem.Load("/home/towa/package/allpix/install/lib/libAllpixObjects.so")

# --- データクラス定義 ---
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

# ★ 新規追加: PropagatedCharge用のデータクラス
@dataclass(slots=True)
class PropagatedCharge:
    charge_type: str  # 'electron' or 'hole'
    global_time: float
    mc_particle: MCParticle
    local_creation_pos: np.ndarray

@dataclass(slots=True)
class Cluster:
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
    def get_position(self, model) -> np.ndarray:
        total_charge = self.charge
        if total_charge == 0:
            return model.get_pixel_center(self.seed_pixel_hit.pixel_index_x, self.seed_pixel_hit.pixel_index_y)
        pos_x, pos_y, pos_z = 0.0, 0.0, 0.0
        for hit in self.pixel_hits:
            center = model.get_pixel_center(hit.pixel_index_x, hit.pixel_index_y)
            pos_x += center[0] * hit.signal
            pos_y += center[1] * hit.signal
            pos_z += center[2] * hit.signal
        return np.array([pos_x / total_charge, pos_y / total_charge, pos_z / total_charge])
    def get_mc_particles(self) -> List[MCParticle]:
        all_particles = [p for hit in self.pixel_hits for p in hit.mc_particles]
        unique_particles = {p.particle_id: p for p in all_particles}
        return list(unique_particles.values())

class DetectorModel:
    def __init__(self, pixel_size_x_um, pixel_size_y_um):
        self.pixel_size = np.array([pixel_size_x_um, pixel_size_y_um])
    def get_pixel_center(self, ix, iy):
        return np.array([(ix + 0.5) * self.pixel_size[0], (iy + 0.5) * self.pixel_size[1], 0.0])
    def get_pixel_index(self, position):
        ix = math.floor(position[0] / self.pixel_size[0])
        iy = math.floor(position[1] / self.pixel_size[1])
        return ix, iy

class AnalysisPixelModule:
    FOUR_NEIGHBORS = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    def __init__(self, config):
        self.config = config
        
        pitch_x = config.get("pixel_pitch_x", 15)
        pitch_y = config.get("pixel_pitch_y", 15)
        self.detector_name = config.get("detector_name", "CE65")
        self.detector_model = DetectorModel(pixel_size_x_um=pitch_x, pixel_size_y_um=pitch_y)
        self.seed_threshold = config.get("seed_threshold", 1000)
        self.neighbor_threshold = config.get("neighbor_threshold", 500)
        self.histograms = {}

        self.fill_3d = config.get("fill_3d", False)

        self.input_file = ROOT.TFile.Open(config["file_name"])
        self.pixel_tree = self.input_file.Get("PixelHit")
        self.mcp_tree = self.input_file.Get("MCParticle")
        self.propagated_tree = self.input_file.Get("PropagatedCharge")
        
        if not self.pixel_tree or not self.mcp_tree or not self.propagated_tree:
             raise RuntimeError("Required TTrees (PixelHit, MCParticle, PropagatedCharge) not found.")

        self.n_entries = self.pixel_tree.GetEntries()
        #self.n_entries = 2000
        
        # sycronize all tree
        self.pixel_tree.AddFriend(self.mcp_tree)
        self.pixel_tree.AddFriend(self.propagated_tree)

        self.counter_names = [ "Total Events", "Skipped: No pixel hits", "Skipped: Multiple primary particles", "Clusters Checked", "Skipped: Hit at EDGE events", "Skipped: Cluster has no primary particle", "Skipped: Residual > 40 um", "Clusters Accepted" ]
        
        print(f"File '{config['file_name']}' opened with {self.n_entries} events.")

    def setup_histograms(self):
        print("Creating histograms...")
        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        inpixel_bins_x = 50
        inpixel_bins_y = 50
        max_cluster_charge_ke = self.config.get("max_cluster_charge_ke", 50.0)

        # (既存のヒストグラム定義は変更なし)
        #self.histograms["inPixel_cluster_size"] = ROOT.TProfile2D("inPixel_cluster_size", ";x/pitch [um];y/pitch [um];cluster size", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        #self.histograms["cluster_charge"] = ROOT.TH1D("cluster_charge", ";charge [ke];counts", 100, 0, max_cluster_charge_ke)
        self.histograms["inPixel_electron_driftTime_90p"] = ROOT.TProfile2D("inPixel_electron_driftTime_90p", ";x/pitch [um];y/pitch [um];90% electron drift time [ns]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["electron_driftTime_90p"] = ROOT.TH1D("electron_driftTime_90p", ";90% electron drift time [ns];counts", 20000, 0, 2)
        self.histograms["inPixel_electron_driftTime_90p_1"] = ROOT.TProfile2D("inPixel_electron_driftTime_90p_1", ";x/pitch [um];y/pitch [um];90% electron drift time [ns]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["electron_driftTime_90p_1"] = ROOT.TH1D("electron_driftTime_90p_1", ";90% electron drift time [ns];counts", 20000, 0, 2)
        
        self.histograms["drift_time_vs_distance"] = ROOT.TH2D(
            "drift_time_vs_distance",
            ";distance to pixel center [um];drift time [ns];counts",
            200, 0, 100,
            1000, 0, 0.5
        )
        self.histograms["cluster_charge_vs_drift_time"] = ROOT.TH2D(
            "cluster_charge_vs_drift_time",
            ";Cluster Charge [ke];90% electron drift time [ns];Counts",
            100, 0, max_cluster_charge_ke,  # Cluster Charge (0-60 ke)
            1000, 0, 0.5                      # Drift Time (0-2 ns)
        )
        
        self.histograms["cluster_size_vs_drift_time"] = ROOT.TH2D(
            "cluster_size_vs_drift_time",
            ";Cluster Size [pixels];90% electron drift time [ns];Counts",
            20, 0, 20,   # Cluster Size (0-20 pixels)
            1000, 0, 0.5   # Drift Time (0-2 ns)
        )

        if self.fill_3d:
            half_pitch_x = self.detector_model.pixel_size[0] / 2.0
            half_pitch_y = self.detector_model.pixel_size[1] / 2.0
            sensor_half_thickness = 25.0

            self.histograms["drift_time_map_xyz"] = ROOT.TProfile3D(
                "drift_time_map_xyz",
                ";In-Pixel X [um];In-Pixel Y [um];Depth Z [um];Average Drift Time [ns]",
                20, -half_pitch_x, half_pitch_x,  # Xビニング (ピクセル内)
                20, -half_pitch_y, half_pitch_y,  # Yビニング (ピクセル内)
                50, -sensor_half_thickness, sensor_half_thickness # Zビニング (深さ)
            )

        n_counters = len(self.counter_names)
        self.histograms["counters"] = ROOT.TH1D("counters", "Event Summary;category;counts", n_counters, 0, n_counters)
        for i, name in enumerate(self.counter_names):
            self.histograms["counters"].GetXaxis().SetBinLabel(i+1, name)

    def do_clustering(self, pixel_hits: List[PixelHit]) -> List[Cluster]:
        # (この関数は変更ありません)
        clusters = []
        hit_map = {(p.pixel_index_x, p.pixel_index_y): p for p in pixel_hits}
        used_pixels = set()
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

    def _fill_histograms_from_buffer(self, buffer: Dict[str, list]):
        # ★ 全ヒストグラムを動的に処理するように改善
        for key, histo in self.histograms.items():
            if not buffer.get(key): continue
            
            if isinstance(histo, ROOT.TProfile2D):
                for x, y, val in buffer[key]:
                    histo.Fill(x, y, val)
            elif isinstance(histo, ROOT.TH1D):
                for val in buffer[key]:
                    histo.Fill(val)
            elif isinstance(histo, ROOT.TH2D):
                for x, y in buffer[key]:
                    histo.Fill(x, y)
            elif isinstance(histo, ROOT.TProfile3D):
                for x, y, z, val in buffer[key]:
                    histo.Fill(x, y, z, val)
        
        for key in buffer:
            buffer[key].clear()

    def run_analysis(self):
        self.setup_histograms()
        branch_name = self.detector_name
        
        N_PIXELS_X, N_PIXELS_Y, MARGIN = 48, 24, 2
        
        self.counters = {name: 0 for name in self.counter_names}
        self.counters["Total Events"] = self.n_entries

        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        offset_vec = np.array([pitch_x / 2, pitch_y / 2, 0])

        BATCH_SIZE = 1000
        buffer = {hist_name: [] for hist_name in self.histograms}

        h_drift_time_vs_dist = self.histograms.get("drift_time_vs_distance")
        h_drift_time_map_xyz = self.histograms.get("drift_time_map_xyz") if self.fill_3d else None

        #for i in range(self.n_entries), desc="Processing Events"):
        for i in range(self.n_entries):
            self.pixel_tree.GetEntry(i)
            
            # --- 1. データデコード ---
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
                mcp.parent = mc_particles.get(mcp.parent)
            
            pixel_objects = getattr(self.pixel_tree, branch_name)
            if not pixel_objects:
                self.counters["Skipped: No pixel hits"] += 1
                continue
            
            pixel_hits = [
                PixelHit(
                    signal=hit_obj.getSignal(),
                    pixel_index_x=hit_obj.getPixel().getIndex().X(),
                    pixel_index_y=hit_obj.getPixel().getIndex().Y(),
                    mc_particles=[mc_particles[p.GetUniqueID()] for p in hit_obj.getMCParticles() if p.GetUniqueID() in mc_particles]
                ) for hit_obj in pixel_objects
            ]

            propagated_objects = getattr(self.propagated_tree, branch_name, [])
            charge_map: Dict[int, List[PropagatedCharge]] = {}
            for prop_obj in propagated_objects:
                mc_particle_ptr = prop_obj.getMCParticle()
                if not mc_particle_ptr: continue
                
                parent_mcp = mc_particles.get(mc_particle_ptr.GetUniqueID())
                if not parent_mcp: continue
                
                #charge_type = 'electron' if prop_obj.getType() == ROOT.allpix.CarrierType.ELECTRON else 'hole'
                if prop_obj.getType() == 255: #ROOT.allpix.CarrierType.ELECTRON = -1
                    charge_type = 'electron'
                elif prop_obj.getType() == 100: #ROOT.allpix.CarrierType.HOLE = 1
                    charge_type = 'hole'
                else:
                    # print(f"prop_obj is {prop_obj.getType()}")
                    # print(f"electron ? {ROOT.allpix.CarrierType.ELECTRON}")
                    # print(f"hole ? {ROOT.allpix.CarrierType.HOLE}")
                    print("There are unkown object number")
                    continue

                pos = prop_obj.getLocalPosition()
                creation_pos_um = np.array([pos.X() * 1000, pos.Y() * 1000, pos.Z() * 1000])

                prop_charge = PropagatedCharge(
                    charge_type=charge_type,
                    global_time=prop_obj.getGlobalTime(),
                    mc_particle=parent_mcp,
                    local_creation_pos=creation_pos_um
                )
                charge_map.setdefault(parent_mcp.particle_id, []).append(prop_charge)

            # --- 2. 解析ロジック ---
            primary_particles = self.get_primary_particles(mc_particles)
            if len(primary_particles) > 1:
                self.counters["Skipped: Multiple primary particles"] += 1
                continue

            clusters = self.do_clustering(pixel_hits)
            for clus in clusters:
                self.counters["Clusters Checked"] += 1
                
                # (既存のクラスタ解析ロジックは変更なし)
                seed_x, seed_y = clus.seed_pixel_hit.pixel_index_x, clus.seed_pixel_hit.pixel_index_y
                if not (MARGIN <= seed_x < N_PIXELS_X - MARGIN and MARGIN <= seed_y < N_PIXELS_Y - MARGIN):
                    self.counters["Skipped: Hit at EDGE events"] += 1
                    continue
                
                cluster_particle_ids = {p.particle_id for p in clus.get_mc_particles()}
                intersection = [p for p in primary_particles if p.particle_id in cluster_particle_ids]
                if not intersection:
                    self.counters["Skipped: Cluster has no primary particle"] += 1
                    continue

                particle = intersection[0]
                particle_pos = particle.local_reference_point + offset_vec
                # cluster_pos = clus.get_position(self.detector_model)
                # residual_r = np.linalg.norm((particle_pos - cluster_pos)[:2])

                # if residual_r > 40:
                #     self.counters["Skipped: Residual > 40 um"] += 1
                #     continue
                
                self.counters["Clusters Accepted"] += 1
                ix, iy = self.detector_model.get_pixel_index(particle_pos)
                pixel_center = self.detector_model.get_pixel_center(ix, iy)
                in_pixel_pos = particle_pos - pixel_center

                for hit in clus.pixel_hits:
                    pixel_center_pos = self.detector_model.get_pixel_center(
                        hit.pixel_index_x,
                        hit.pixel_index_y
                    )

                    for mcp in hit.mc_particles:
                        propagated_charges = charge_map.get(mcp.particle_id, [])

                        for pc in propagated_charges:
                            if pc.charge_type == 'electron':
                                drift_time_ns = pc.global_time * 1e-3
                                creation_pos = pc.local_creation_pos
                                distance_um = np.linalg.norm(creation_pos - pixel_center_pos)
                                
                                if h_drift_time_vs_dist:
                                    h_drift_time_vs_dist.Fill(distance_um, drift_time_ns)

                                #buffer["drift_time_vs_distance"].append((distance_um, drift_time_ns))

                                if self.fill_3d:
                                    relative_pos = creation_pos - pixel_center_pos 
                                    # buffer["drift_time_map_xyz"].append((
                                    #     relative_pos[0], # ピクセル内X座標
                                    #     relative_pos[1], # ピクセル内Y座標
                                    #     relative_pos[2], # センサ内Z座標 (深さ)
                                    #     drift_time_ns    # ドリフト時間
                                    # ))
                                    h_drift_time_map_xyz.Fill(
                                        relative_pos[0],
                                        relative_pos[1],
                                        relative_pos[2],
                                        drift_time_ns
                                    )
                # propagated_charges_for_particle = charge_map.get(particle.particle_id, [])
                # if propagated_charges_for_particle:
                #     electron_drift_times_ps = [pc.global_time for pc in propagated_charges_for_particle if pc.charge_type == 'electron']
                    
                #     if electron_drift_times_ps:
                #         # 90パーセンタイル値を計算 (ps -> nsに変換)
                #         time_90_percent_ns = np.percentile(electron_drift_times_ps, 90) * 1e-3

                #         #print(time_90_percent_ns)
                        
                #         buffer["electron_driftTime_90p_1"].append(time_90_percent_ns)
                #         #print(time_90_percent_ns)
                #         buffer["inPixel_electron_driftTime_90p_1"].append((in_pixel_pos[0], in_pixel_pos[1], time_90_percent_ns))

                electron_drift_time_ps = [
                    pc.global_time
                    for hit in clus.pixel_hits
                    for mcp in hit.mc_particles
                    for pc in charge_map.get(mcp.particle_id, [])
                    if pc.charge_type == 'electron'
                ]

                if electron_drift_time_ps:
                    time_90_percent_ns = np.percentile(electron_drift_time_ps, 90) * 1e-3
                    buffer["electron_driftTime_90p"].append(time_90_percent_ns)
                    buffer["inPixel_electron_driftTime_90p"].append((in_pixel_pos[0], in_pixel_pos[1], time_90_percent_ns))
                    buffer["cluster_charge_vs_drift_time"].append((clus.charge / 1000.0, time_90_percent_ns)) # e -> ke に変換
                    buffer["cluster_size_vs_drift_time"].append((clus.size, time_90_percent_ns))

            if (i + 1) % BATCH_SIZE == 0:
                self._fill_histograms_from_buffer(buffer)
        
        if any(buffer.values()):
            self._fill_histograms_from_buffer(buffer)

        print("Analysis finished.")
        # (レポート表示は変更なし)
        print("----------------- Report -----------------")
        for key in self.counter_names:
            print(f"{key:<40}:{self.counters.get(key, 0)}")
        print("------------------------------------------")

    def finalize(self):
        for i, name in enumerate(self.counter_names):
            count = self.counters.get(name, 0)
            self.histograms["counters"].SetBinContent(i + 1, count)
        
        output_filename = self.config.get("output_file_name", "analysis_py.root")
        output_file = ROOT.TFile(output_filename, "RECREATE")
        for histo in self.histograms.values():
            histo.Write()
        output_file.Close()
        print(f"Histograms written to {output_filename}")


if __name__ == '__main__':
    # (mainブロックは変更なし)
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Allpix Squared Analysis script.")
    # ... args ...
    parser.add_argument("-i", "--input_file", type=str, help="Path to the input ROOT file from Allpix Squared Simulation Data.")
    parser.add_argument("-o", "--output", default="analysis_py.root", type=str, help="Path for the output ROOT file.")
    parser.add_argument("-st", "--seed_threshold", type=int, default=0, help="Seed_threshold for analysis (electron num).")
    parser.add_argument("-nt", "--neighbor_threshold", type=int, default=0, help="Neighbor_threshold for analysis (electron num).")
    parser.add_argument("-p", "--pixel_pitch", type=float, default=22.5, help="pixel pitch about DUT")
    parser.add_argument("-v", "--voltage", type=int, default=10, help="Chip voltage about DUT (e.g., 10, 7, 4)")
    parser.add_argument("-b", "--beam_type", type=str, default="e3GeV", help="Beam information e.g., e3GeV")
    parser.add_argument("-m", "--model", type=str, default="masetti", help="Model name for Electron calculation")
    parser.add_argument("-n", "--name", type=str, default="CE65", help="DUT name")
    parser.add_argument(
        "--fill3D", 
        action="store_true", 
        help="Enable filling of the 3D drift time map (can be slow)."
    )

    args = parser.parse_args()    
    config = {
        "file_name": args.input_file,
        "output_file_name": args.output,
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
        "max_cluster_charge_ke": 60.0,
        "fill_3d": args.fill3D
    }

    try:
        analyzer = AnalysisPixelModule(config)
        analyzer.run_analysis()
        analyzer.finalize()
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

    end_time = time.time()
    # ... (時間表示は変更なし) ...
    print(f"\nTotal time: {end_time - start_time:.2f} seconds")