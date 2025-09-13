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

gSystem.Load("/home/towa/package/allpix/install/lib/libAllpixObjects.so")

# C++のPixelHit, MCParticleなどのデータ構造を模倣するためのデータクラス
# (この部分は変更ありません)
@dataclass
class PixelHit:
    signal: float
    pixel_index_x: int
    pixel_index_y: int
    mc_particles: List[Any]

@dataclass
class MCParticle:
    particle_id: int
    parent: Any
    local_reference_point: np.ndarray

@dataclass
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
    
    def get_position(self, model) -> np.ndarray:
        total_charge = self.charge
        if total_charge == 0:
            return model.get_pixel_center(
                self.seed_pixel_hit.pixel_index_x,
                self.seed_pixel_hit.pixel_index_y
            )
        
        pos_x, pos_y, pos_z = 0.0, 0.0, 0.0
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

        # PyROOTでファイルとTTreeを開く
        self.input_file = ROOT.TFile.Open(config["file_name"])
        self.pixel_tree = self.input_file.Get("PixelHit")
        self.mcp_tree = self.input_file.Get("MCParticle")
        self.n_entries = self.pixel_tree.GetEntries()
        #self.n_entries = 100
        
        # 2つのTTreeが同じイベントを指すように同期させる
        self.pixel_tree.AddFriend(self.mcp_tree)

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
        
        print(f"File '{config['file_name']}' opened with {self.n_entries} events.")

    def setup_histograms(self):
        # (この関数は変更ありません)
        print("Creating histograms...")
        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        inpixel_bins_x = 50
        inpixel_bins_y = 50
        max_cluster_charge_ke = self.config.get("max_cluster_charge_ke", 50.0)
        self.histograms["inPixel_cluster_size"] = ROOT.TProfile2D("inPixel_cluster_size", ";x/pitch [um];y/pitch [um];cluster size", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_cluster_charge"] = ROOT.TProfile2D("inPixel_cluster_charge", ";x/pitch [um];y/pitch [um];cluster charge [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_seed_charge"] = ROOT.TProfile2D("inPixel_seed_charge", ";x/pitch [um];y/pitch [um];seed charge [ke]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_r"] = ROOT.TProfile2D("inPixel_residual_r", ";x/pitch [um];y/pitch [um];residual r [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_x"] = ROOT.TProfile2D("inPixel_residual_x", ";x/pitch [um];y/pitch [um];residual x [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_y"] = ROOT.TProfile2D("inPixel_residual_y", ";x/pitch [um];y/pitch [um];residual y [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        self.histograms["inPixel_residual_xy2"] = ROOT.TProfile2D("inPixel_residual_xy2", ";x/pitch [um];y/pitch [um];x + y / 2 [um]", inpixel_bins_x, -pitch_x / 2, pitch_x / 2, inpixel_bins_y, -pitch_y / 2, pitch_y / 2)
        # self.histograms["inPixel_residual"] = ROOT.TProfile2D("inPixel_residual", ";x/pitch [um];y/pitch [um];residual [um]", inpixel_bins_x, -1000, 1000, inpixel_bins_y, -1000, 1000)
        #self.histograms["test"] = ROOT.TH2D("test", "test", 1000, -50, 50, 1000, -50, 50)
        self.histograms["cluster_charge"] = ROOT.TH1D("cluster_charge", ";charge [ke];counts", 1000, 0, max_cluster_charge_ke)
        self.histograms["cluster_size"] = ROOT.TH1D("cluster_size", ";cluster size;counts", 20, 0.5, 20.5)
        self.histograms["seed_charge"] = ROOT.TH1D("seed_charge", ";charge [ke];counts", 1000, 0, max_cluster_charge_ke)
        self.histograms["residual_x"] = ROOT.TH1D("residual_x", ";residual x [um];counts", 10000, -40, 40)
        self.histograms["residual_y"] = ROOT.TH1D("residual_y", ";residual y [um];counts", 10000, -40, 40)
        self.histograms["residual_r"] = ROOT.TH1D("residual_r", ";residual r [um];counts", 10000, 0, 40)

        n_counters = len(self.counter_names)
        self.histograms["counters"] = ROOT.TH1D("counters", "Event Summary;category;counts", n_counters, 0, n_counters)
        for i, name in enumerate(self.counter_names):
            self.histograms["counters"].GetXaxis().SetBinLabel(i+1, name)

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
        for val in buffer["cluster_charge"]:
            self.histograms["cluster_charge"].Fill(val)
        for val in buffer["cluster_size"]:
            self.histograms["cluster_size"].Fill(val)
        for val in buffer["seed_charge"]:
            self.histograms["seed_charge"].Fill(val)
        for val in buffer["residual_x"]:
            self.histograms["residual_x"].Fill(val)
        for val in buffer["residual_y"]:
            self.histograms["residual_y"].Fill(val)
        for val in buffer["residual_r"]:
            self.histograms["residual_r"].Fill(val)

        for x, y, val in buffer["inPixel_cluster_size"]:
            self.histograms["inPixel_cluster_size"].Fill(x, y, val)
        for x, y, val in buffer["inPixel_residual_r"]:
            self.histograms["inPixel_residual_r"].Fill(x, y, val)
        for x, y, val in buffer["inPixel_residual_x"]:
            self.histograms["inPixel_residual_x"].Fill(x, y, val)
        for x, y, val in buffer["inPixel_residual_y"]:
            self.histograms["inPixel_residual_y"].Fill(x, y, val)
        for x, y, val in buffer["inPixel_residual_xy2"]:
            self.histograms["inPixel_residual_xy2"].Fill(x, y, val)
        for x, y, val in buffer["inPixel_seed_charge"]:
            self.histograms["inPixel_seed_charge"].Fill(x, y, val)
        for x, y, val in buffer["inPixel_cluster_charge"]:
            self.histograms["inPixel_cluster_charge"].Fill(x, y, val)
            
        for key in buffer:
            buffer[key].clear()

    def run_analysis(self):
        """PyROOTを使ってイベントループを実行し、解析を行う"""
        self.setup_histograms()
        
        branch_name = "CE65"

        N_PIXELS_X = 48
        N_PIXELS_Y = 24
        MARGIN = 2

        # counters = {
        #     "Total Events": self.n_entries,
        #     "Skipped: No pixel hits": 0,
        #     "Skipped: Multiple primary particles": 0,
        #     "Clusters Checked": 0,
        #     "Skipped: Cluster has no primary particle": 0,
        #     "Skipped: Residual > 40 um": 0,
        #     "Clusters Accepted": 0,
        # }
        self.counters = {name: 0 for name in self.counter_names}
        self.counters["Total Events"] = self.n_entries

        pitch_x = self.detector_model.pixel_size[0]
        pitch_y = self.detector_model.pixel_size[1]
        offset_x = pitch_x / 2
        offset_y = pitch_y / 2

        BATCH_SIZE =1000
        buffer = {hist_name: [] for hist_name in self.histograms}

        #for i in tqdm(range(self.n_entries)):
        for i in range(self.n_entries):
            self.pixel_tree.GetEntry(i)
            
            # --- 1. データデコード ---
            
            # PyROOTはC++のstd::vectorを返す。これは直接ループ処理できる。
            mcp_objects = getattr(self.mcp_tree, branch_name)
            
            mc_particles: Dict[int, MCParticle] = {}
            # まずPythonのMCParticleオブジェクトを全作成
            for mcp_obj in mcp_objects:
                pid = mcp_obj.GetUniqueID()
                parent_ptr = mcp_obj.getParent()
                parent_id = parent_ptr.GetUniqueID() if parent_ptr else 0
                
                ref_point = mcp_obj.getLocalReferencePoint()
                mc_particles[pid] = MCParticle(
                    particle_id=pid,
                    parent=parent_id, # まずはIDを入れておく
                    local_reference_point=np.array([ref_point.X()*1000, ref_point.Y()*1000, ref_point.Z()*1000])
                )

            # 親への参照を解決
            for mcp in mc_particles.values():
                parent_id = mcp.parent
                if parent_id != 0 and parent_id in mc_particles:
                    mcp.parent = mc_particles[parent_id]
                else:
                    mcp.parent = None
            
            # PixelHitオブジェクトを取得
            pixel_objects = getattr(self.pixel_tree, branch_name)
            if len(pixel_objects) == 0:
                self.counters["Skipped: No pixel hits"] += 1
                continue
            
            pixel_hits = []
            if len(pixel_objects) == 0:
                continue

            for hit_obj in pixel_objects:
                # このヒットに関連するMCParticleをPython辞書から見つける
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

            # --- 2. 解析ロジック (ここからは以前のコードとほぼ同じ) ---
            primary_particles = self.get_primary_particles(mc_particles)
            if len(primary_particles) > 1:
                self.counters["Skipped: Multiple primary particles"] += 1
                continue

            clusters = self.do_clustering(pixel_hits)

            for clus in clusters:
                self.counters["Clusters Checked"] += 1
                # self.histograms["cluster_charge"].Fill(clus.charge / 1000.0)
                # self.histograms["cluster_size"].Fill(clus.size)
                # self.histograms["seed_charge"].Fill(clus.seed_pixel_hit.signal / 1000.0)

                seed_x = clus.seed_pixel_hit.pixel_index_x
                seed_y = clus.seed_pixel_hit.pixel_index_y

                if not (MARGIN <= seed_x < N_PIXELS_X - MARGIN and MARGIN <= seed_y < N_PIXELS_Y - MARGIN):
                    self.counters["Skipped: Hit at EDGE events"] += 1
                    continue
                
                #cluster_mc_particles = clus.get_mc_particles()
                #intersection = [p for p in primary_particles if p in cluster_mc_particles]
                cluster_particles_ids = {p.particle_id for p in clus.get_mc_particles()}
                #cluster_particles_set = set(clus.get_mc_particles())
                intersection = [p for p in primary_particles if p.particle_id in cluster_particles_ids]
                
                if not intersection:
                    self.counters["Skipped: Cluster has no primary particle"] += 1
                    continue

                particle = intersection[0]
                particle_pos = particle.local_reference_point + np.array([offset_x, offset_y, 0])
                cluster_pos = clus.get_position(self.detector_model)

                # ix, iy = self.detector_model.get_pixel_index(particle_pos)
                # pixel_center = self.detector_model.get_pixel_center(ix, iy)
                # in_pixel_pos = particle_pos - pixel_center
                
                residual_vec = particle_pos - cluster_pos
                residual_r = np.linalg.norm(residual_vec[:2])

                # print("----------- MCParticle Hit Position -----------")
                # print(particle_pos)
                # print("----------- Clustering Hit Position -----------")
                # print(cluster_pos)
                # print("----------- Residual Value -----------")
                # print(residual_vec)
                # print()
                # print("==================================================")

                if residual_r > 40:
                    self.counters["Skipped: Residual > 40 um"] += 1
                    continue

                self.counters["Clusters Accepted"] += 1

                # self.histograms["cluster_charge"].Fill(clus.charge / 1000.0)
                # self.histograms["cluster_size"].Fill(clus.size)
                # self.histograms["seed_charge"].Fill(clus.seed_pixel_hit.signal / 1000.0)

                ix, iy = self.detector_model.get_pixel_index(particle_pos)
                pixel_center = self.detector_model.get_pixel_center(ix, iy)
                in_pixel_pos = particle_pos - pixel_center
                #print(in_pixel_pos)

                # ix_reco, iy_reco = self.detector_model.get_pixel_index(cluster_pos)

                # reference_pixel_center = self.detector_model.get_pixel_center(
                #     ix_reco,
                #     iy_reco
                # )

                # relative_pos = particle_pos - reference_pixel_center

                # in_pixel_pos_x = (relative_pos[0] + pitch_x / 2) % pitch_x - pitch_x / 2
                # in_pixel_pos_y = (relative_pos[1] + pitch_y / 2) % pitch_y - pitch_y / 2

                # self.histograms["cluster_charge"].Fill(clus.charge / 1000.0)
                # self.histograms["cluster_size"].Fill(clus.size)
                # self.histograms["seed_charge"].Fill(clus.seed_pixel_hit.signal / 1000.0)
                # self.histograms["residual_x"].Fill(residual_vec[0])
                # self.histograms["residual_y"].Fill(residual_vec[1])
                # self.histograms["residual_r"].Fill(residual_r)
                # self.histograms["inPixel_cluster_size"].Fill(in_pixel_pos[0], in_pixel_pos[1], clus.size)
                # self.histograms["inPixel_residual_r"].Fill(in_pixel_pos[0], in_pixel_pos[1], residual_r)
                # self.histograms["inPixel_residual_x"].Fill(in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[0]))
                # self.histograms["inPixel_residual_y"].Fill(in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[1]))
                # self.histograms["inPixel_residual_xy2"].Fill(in_pixel_pos[0], in_pixel_pos[1], (abs(residual_vec[0]) + abs(residual_vec[1])) / 2)
                # self.histograms["inPixel_seed_charge"].Fill(in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal / 1000.0)
                # self.histograms["inPixel_cluster_charge"].Fill(in_pixel_pos[0], in_pixel_pos[1], clus.charge / 1000.0)
                #self.histograms["test"].Fill(in_pixel_pos[0], in_pixel_pos[1])

                # pitch_x = self.detector_model.pixel_size[0]
                # pitch_y = self.detector_model.pixel_size[1]

                # print(f"\n--- DEBUG (Event {i}, Cluster Accepted) ---")
                # print(f"  Particle Pos (x, y):  ({particle_pos[0]:.3f}, {particle_pos[1]:.3f})")
                
                # ix, iy = self.detector_model.get_pixel_index(particle_pos)
                # print(f"  Pixel Index (ix, iy): ({ix}, {iy})")
                
                # pixel_center = self.detector_model.get_pixel_center(ix, iy)
                # print(f"  Pixel Center (x, y):  ({pixel_center[0]:.3f}, {pixel_center[1]:.3f})")

                # in_pixel_pos = particle_pos - pixel_center
                # print(f"  => InPixel Pos (x, y): ({in_pixel_pos[0]:.3f}, {in_pixel_pos[1]:.3f})")
                # print(f"  Histogram Range X:    ({-pitch_x / 2:.3f} to {pitch_x / 2:.3f})")
                # print(f"  Histogram Range Y:    ({-pitch_y / 2:.3f} to {pitch_y / 2:.3f})")
                # print(f"------------------------------------")

                buffer["cluster_charge"].append(clus.charge / 1000.0)
                buffer["cluster_size"].append(clus.size)
                buffer["seed_charge"].append(clus.seed_pixel_hit.signal / 1000.0)
                buffer["residual_x"].append(residual_vec[0])
                buffer["residual_y"].append(residual_vec[1])
                buffer["residual_r"].append(residual_r)
                buffer["inPixel_cluster_size"].append((in_pixel_pos[0], in_pixel_pos[1], clus.size))
                buffer["inPixel_residual_r"].append((in_pixel_pos[0], in_pixel_pos[1], residual_r))
                buffer["inPixel_residual_x"].append((in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[0])))
                buffer["inPixel_residual_y"].append((in_pixel_pos[0], in_pixel_pos[1], abs(residual_vec[1])))
                buffer["inPixel_residual_xy2"].append((in_pixel_pos[0], in_pixel_pos[1], (abs(residual_vec[0]) + abs(residual_vec[1])) / 2))
                buffer["inPixel_seed_charge"].append((in_pixel_pos[0], in_pixel_pos[1], clus.seed_pixel_hit.signal / 1000.0))
                buffer["inPixel_cluster_charge"].append((in_pixel_pos[0], in_pixel_pos[1], clus.charge / 1000.0))
                # --------------------------------------------------------

            if (i + 1) % BATCH_SIZE == 0:
                self._fill_histograms_from_buffer(buffer)
            
        if any(buffer.values()):
            self._fill_histograms_from_buffer(buffer)

        print("Analysis finished.")
        print("----------------- Report -----------------")
        for key, value in self.counters.items():
            print(f"{key:<40}:{value}")
        print("------------------------------------------")
        

    def finalize(self):
        for i, name in enumerate(self.counter_names):
            count = self.counters.get(name, 0)
            self.histograms["counters"].SetBinContent(i + 1, count)
        
        output_filename = self.config.get("output_file_name", "analysis_py.root")
        output_file = ROOT.TFile(output_filename, "RECREATE")
        for name, histo in self.histograms.items():
            histo.Write()
        output_file.Close()
        print(f"Histograms written to {output_filename}")


if __name__ == '__main__':
    start_time = time.time()
    
    parser = argparse.ArgumentParser(
        description="Allpix Squared Analysis script.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i", "--input_file", type=str, help="Path to the input ROOT file from Allpix Squared Simulation Data.")
    parser.add_argument("-o", "--output", default="analysis_py.root", type=str, help="Path for the output ROOT file.")
    parser.add_argument("-st", "--seed_threshold", type=int, default=0, help="Seed_threshold for analysis (electron num).")
    parser.add_argument("-nt", "--neighbor_threshold", type=int, default=0, help="Neighbor_threshold for analysis (electron num).")
    parser.add_argument("-p", "--pixel_pitch", type=float, default=22.5, help="pixel pitch about DUT")
    parser.add_argument("-v", "--voltage", type=int, default=10, help="Chip voltage about DUT (e.g., 10, 7, 4)")
    parser.add_argument("-b", "--beam_type", type=str, default="e3GeV", help="Beam information e.g., e3GeV")
    parser.add_argument("-m", "--model", type=str, default="masetti", help="Model name for Electron calculation")
    parser.add_argument("-n", "--name", type=str, default="CE65", help="DUT name")

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
        "max_cluster_charge_ke": 60.0
    }

    try:
        f = ROOT.TFile.Open(config["file_name"])
        if not f or f.IsZombie():
            raise FileNotFoundError
        if not f.Get("PixelHit") or not f.Get("MCParticle"):
             print(f"Error: Input file '{config['file_name']}' is missing required TTrees.")
             print(f"Required: 'PixelHit' and 'MCParticle'")
             exit()
        f.Close()
        print(f"Using existing file: {config['file_name']}")
    except (FileNotFoundError, OSError):
        print(f"Error: Input file '{config['file_name']}' not found.")
        exit()

    analyzer = AnalysisPixelModule(config)
    analyzer.run_analysis()
    analyzer.finalize()

    end_time = time.time()
    elapsed_time = end_time - start_time
    start_str = datetime.fromtimestamp(start_time).strftime("%Y/%m/%d %H:%M:%S")
    end_str = datetime.fromtimestamp(end_time).strftime("%Y/%m/%d %H:%M:%S")
    print("------------------------------------------")
    print(f"Start time : {start_str}")
    print(f"End time   : {end_str}")
    print(f"Total time : {elapsed_time:.2f} seconds")

    minutes = int(elapsed_time // 60)
    seconds = elapsed_time % 60
    if minutes > 0:
        print(f"Which is:     {minutes} minute(s) and {seconds:.2f} second(s)")
    print("------------------------------------------")