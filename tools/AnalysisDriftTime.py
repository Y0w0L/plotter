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
import multiprocessing
import os
import traceback
import gc

# Load Allpix Squared objects library
gSystem.Load("/home/towa/package/allpix/install/lib/libAllpixObjects.so")
ROOT.ROOT.DisableImplicitMT()

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
    
    def get_position(self, model, one_bit: bool) -> np.ndarray:
        centers = np.array([model.get_pixel_center(h.pixel_index_x, h.pixel_index_y) for h in self.pixel_hits])
        if one_bit:
            return np.mean(centers, axis=0)
        else:
            signals = np.array([h.signal for h in self.pixel_hits])
            total_charge = np.sum(signals)
            if total_charge <= 0:
                return centers[0]
            return np.average(centers, axis=0, weights=signals)

    def get_mc_particles(self) -> List[Any]:
        unique_particles = {p.particle_id: p for hit in self.pixel_hits for p in hit.mc_particles}
        return list(unique_particles.values())

class DetectorModel:
    def __init__(self, pixel_size_x_um, pixel_size_y_um):
        self.pixel_size = np.array([pixel_size_x_um, pixel_size_y_um])

    def get_pixel_center(self, ix, iy):
        return np.array([(ix + 0.5) * self.pixel_size[0], (iy + 0.5) * self.pixel_size[1], 0.0])

    def get_pixel_index(self, position_um: np.ndarray):
        indices = np.floor(position_um[:2] / self.pixel_size)
        return int(indices[0]), int(indices[1])

class AnalysisPixelModule:
    def __init__(self, config):
        self.config = config
        pitch = np.array([config.get("pixel_pitch_x", 22.5), config.get("pixel_pitch_y", 22.5)])
        self.detector_model = DetectorModel(pitch[0], pitch[1])
        self.seed_threshold = config.get("seed_threshold", 1000)
        self.neighbor_threshold = config.get("neighbor_threshold", 500)
        self.histograms = {}
        self.event_histograms = {}  
        self.one_bit_processing = config.get("one_bit", False)
        self.fill_3d = config.get("fill_3d", False)
        self.save_event_hists = config.get("save_event_hists", False)
        self.save_per_size = config.get("save_per_size", False) # Flag for PerClusterSize
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
            "Total Events", "Skipped: No pixel hits", "Skipped: Multiple primary particles",
            "Clusters Checked", "Skipped: Hit at EDGE events", 
            "Skipped: Cluster has no primary particle", "Skipped: Residual > 40 um", "Clusters Accepted"
        ]
        self.counters = {name: 0 for name in self.counter_names}

    def setup_histograms(self):
        pitch_x, pitch_y = self.detector_model.pixel_size
        bins = 50
        max_q = self.config.get("max_cluster_charge", 20000)
        q_bins = 100

        # In-Pixel Profiles
        for name, title in [
            ("inPixel_efficiency", "Efficiency"), ("inPixel_cluster_size", "Cluster Size"),
            ("inPixel_cluster_charge", "Cluster Charge [ke]"), ("inPixel_seed_charge", "Seed Charge [ke]"),
            ("inPixel_neighbor_charge_sum", "Neighbor Sum [ke]"), ("inPixel_residual_r", "Residual R [um]"),
            ("inPixel_seed_ratio", "Seed Ratio"), ("inPixel_multi_hit_prob", "Multi-hit Prob"),
        ]:
            self.histograms[name] = ROOT.TProfile2D(name, f";x [um];y [um];{title}", bins, -pitch_x/2, pitch_x/2, bins, -pitch_y/2, pitch_y/2)

        self.drift_fractions = [50, 60, 70, 80, 90]
        for f in self.drift_fractions:
            h_name = f"electron_driftTime_{f}p"
            self.histograms[h_name] = ROOT.TH1D(h_name, f";{f}% Drift Time [ns];Counts", 1000, 0, 0.5)
            p_name = f"inPixel_electron_driftTime_{f}p"
            self.histograms[p_name] = ROOT.TProfile2D(p_name, f";x [um];y [um];{f}% Drift Time [ns]", bins, -pitch_x/2, pitch_x/2, bins, -pitch_y/2, pitch_y/2)

        # 1D Spectra
        self.histograms["cluster_charge"] = ROOT.TH1D("cluster_charge", ";Charge [ke];Counts", q_bins, 0, max_q)
        self.histograms["cluster_size"] = ROOT.TH1D("cluster_size", ";Size;Counts", 20, 0.5, 20.5)
        self.histograms["seed_charge"] = ROOT.TH1D("seed_charge", ";Charge [ke];Counts", q_bins, 0, max_q)
        self.histograms["residual_r"] = ROOT.TH1D("residual_r", ";Residual R [um];Counts", 100, 0, 40)
        self.histograms["residual_x"] = ROOT.TH1D("residual_x", ";Residual X [um];Counts", 100, -40, 40)
        self.histograms["residual_y"] = ROOT.TH1D("residual_y", ";Residual Y [um];Counts", 100, -40, 40)
        self.histograms["drift_time_spectrum"] = ROOT.TH1D("drift_time_spectrum", ";All Drift Time [ns];Counts", 1000, 0, 0.5)

        # 2D Correlation
        self.histograms["cluster_charge_vs_drift_time_90p"] = ROOT.TH2D("cluster_charge_vs_drift_time_90p", ";Charge [ke];90% Drift Time [ns]", q_bins, 0, max_q, 200, 0, 0.2)
        self.histograms["cluster_size_vs_drift_time_90p"] = ROOT.TH2D("cluster_size_vs_drift_time_90p", ";Size;90% Drift Time [ns]", 20, 0.5, 20.5, 200, 0, 0.2)

        # Counters
        self.histograms["counters"] = ROOT.TH1D("counters", ";Category;Counts", len(self.counter_names), 0, len(self.counter_names))
        for i, name in enumerate(self.counter_names): self.histograms["counters"].GetXaxis().SetBinLabel(i+1, name)

        # Per-size histograms (Conditional)
        if self.save_per_size:
            for i in range(1, self.max_cluster_size_hist + 2):
                suffix = f"_size_{i}" if i <= self.max_cluster_size_hist else f"_size_{self.max_cluster_size_hist+1}_plus"
                self.histograms[f"seed_charge{suffix}"] = ROOT.TH1D(f"seed_charge{suffix}", ";Charge [ke];Counts", q_bins, 0, max_q)
                self.histograms[f"residual_r{suffix}"] = ROOT.TH1D(f"residual_r{suffix}", ";R [um];Counts", 100, 0, 40)
                self.histograms[f"residual_x{suffix}"] = ROOT.TH1D(f"residual_x{suffix}", ";X [um];Counts", 100, -40, 40)
                self.histograms[f"residual_y{suffix}"] = ROOT.TH1D(f"residual_y{suffix}", ";Y [um];Counts", 100, -40, 40)

        if self.fill_3d:
            self.histograms["drift_time_map_xyz"] = ROOT.TProfile3D("drift_time_map_xyz", ";X [um];Y [um];Z [um];Time [ns]", 10, -pitch_x/2, pitch_x/2, 10, -pitch_y/2, pitch_y/2, 25, -25, 25)

    def do_clustering(self, hits: List[PixelHit]) -> List[Cluster]:
        clusters = []
        used = set()
        hit_map = {(h.pixel_index_x, h.pixel_index_y): h for h in hits}
        for seed in sorted(hits, key=lambda h: h.signal, reverse=True):
            sid = (seed.pixel_index_x, seed.pixel_index_y)
            if sid in used or seed.signal < self.seed_threshold: continue
            cluster = Cluster(seed)
            used.add(sid)
            queue = [seed]
            while queue:
                curr = queue.pop(0)
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    nid = (curr.pixel_index_x + dx, curr.pixel_index_y + dy)
                    if nid in hit_map and nid not in used and hit_map[nid].signal >= self.neighbor_threshold:
                        cluster.add_pixel_hit(hit_map[nid])
                        used.add(nid)
                        queue.append(hit_map[nid])
            clusters.append(cluster)
        return clusters

    def _fill_histograms_from_buffer(self, buffer: Dict[str, list]):
        for key, data_list in buffer.items():
            if not data_list or key not in self.histograms: continue
            histo = self.histograms[key]
            for data in data_list:
                if isinstance(data, tuple): histo.Fill(*data)
                else: histo.Fill(data)
            data_list.clear()

    def run_analysis(self, events_per_chunk=500):
        self.setup_histograms()
        branch = self.config.get("detector_name", "CE65")
        N_X, N_Y, MARGIN = 48, 24, 2
        batch_size = 100
        buffer = {k: [] for k in self.histograms}
        pitch = self.detector_model.pixel_size
        offset_vec = np.array([pitch[0]/2, pitch[1]/2, 0.0])

        for i_ev in tqdm(range(self.start_entry, self.end_entry), position=self.worker_id, leave=False):
            if self.pixel_tree.GetEntry(i_ev) <= 0: continue
            self.counters["Total Events"] += 1

            mcp_objs = getattr(self.mcp_tree, branch)
            mc_particles, primaries = {}, []
            for obj in mcp_objs:
                pid = obj.GetUniqueID()
                ref = obj.getLocalReferencePoint()
                mcp = MCParticle(pid, obj.getParent().GetUniqueID() if obj.getParent() else 0, np.array([ref.X(), ref.Y(), ref.Z()]) * 1000.0)
                mc_particles[pid] = mcp
                if mcp.parent == 0: primaries.append(mcp)

            if len(primaries) != 1:
                self.counters["Skipped: Multiple primary particles"] += 1
                continue
            
            target_mcp = primaries[0]
            particle_pos = target_mcp.local_reference_point + offset_vec
            ix_mc, iy_mc = self.detector_model.get_pixel_index(particle_pos)
            is_in_roi = (MARGIN <= ix_mc < N_X - MARGIN and MARGIN <= iy_mc < N_Y - MARGIN)
            in_pixel_pos_mc = particle_pos - self.detector_model.get_pixel_center(ix_mc, iy_mc)

            pixel_objs = getattr(self.pixel_tree, branch)
            hits = [PixelHit(o.getSignal(), o.getPixel().getIndex().X(), o.getPixel().getIndex().Y(), 
                             [mc_particles[p.GetUniqueID()] for p in o.getMCParticles() if p.GetUniqueID() in mc_particles]) for o in pixel_objs]
            clusters = self.do_clustering(hits)
            best_cl = max(clusters, key=lambda c: c.seed_pixel_hit.signal) if clusters else None
            
            is_rec, res_vec = False, None
            if best_cl:
                if target_mcp.particle_id in {p.particle_id for p in best_cl.get_mc_particles()}:
                    res_vec = particle_pos - best_cl.get_position(self.detector_model, self.one_bit_processing)
                    if np.linalg.norm(res_vec[:2]) <= 40: is_rec = True

            if is_in_roi:
                buffer["inPixel_efficiency"].append((in_pixel_pos_mc[0], in_pixel_pos_mc[1], 1.0 if is_rec else 0.0))

            if not is_rec:
                if not hits: self.counters["Skipped: No pixel hits"] += 1
                elif best_cl is None: pass
                elif res_vec is None: self.counters["Skipped: Cluster has no primary particle"] += 1
                else: self.counters["Skipped: Residual > 40 um"] += 1
                continue

            seed = best_cl.seed_pixel_hit
            if not (MARGIN <= seed.pixel_index_x < N_X - MARGIN and MARGIN <= seed.pixel_index_y < N_Y - MARGIN):
                self.counters["Skipped: Hit at EDGE events"] += 1
                continue

            self.counters["Clusters Accepted"] += 1
            pix_center = self.detector_model.get_pixel_center(seed.pixel_index_x, seed.pixel_index_y)
            in_x, in_y = (particle_pos - pix_center)[0], (particle_pos - pix_center)[1]
            q_sum, q_seed, cl_size = best_cl.charge, seed.signal, best_cl.size
            res_r = np.linalg.norm(res_vec[:2])

            buffer["cluster_charge"].append(q_sum)
            buffer["cluster_size"].append(cl_size)
            buffer["seed_charge"].append(q_seed)
            buffer["residual_r"].append(res_r)
            buffer["residual_x"].append(res_vec[0])
            buffer["residual_y"].append(res_vec[1])
            buffer["inPixel_cluster_size"].append((in_x, in_y, cl_size))
            buffer["inPixel_seed_charge"].append((in_x, in_y, q_seed))
            buffer["inPixel_cluster_charge"].append((in_x, in_y, q_sum))
            buffer["inPixel_neighbor_charge_sum"].append((in_x, in_y, q_sum - q_seed))
            buffer["inPixel_residual_r"].append((in_x, in_y, res_r))
            if q_sum > 0: buffer["inPixel_seed_ratio"].append((in_x, in_y, q_seed / q_sum))
            buffer["inPixel_multi_hit_prob"].append((in_x, in_y, 1.0 if cl_size > 1 else 0.0))

            # Buffer filling for PerSize (Conditional)
            if self.save_per_size:
                s_suffix = f"_size_{cl_size}" if cl_size <= self.max_cluster_size_hist else f"_size_{self.max_cluster_size_hist+1}_plus"
                if f"seed_charge{s_suffix}" in buffer:
                    buffer[f"seed_charge{s_suffix}"].append(q_seed)
                    buffer[f"residual_r{s_suffix}"].append(res_r)
                    buffer[f"residual_x{s_suffix}"].append(res_vec[0])
                    buffer[f"residual_y{s_suffix}"].append(res_vec[1])

            prop_objs = getattr(self.propagated_tree, branch)
            dts = []
            hist_name = f"event_{i_ev}_drift_time"
            if self.save_event_hists:
                self.event_histograms[hist_name] = ROOT.TH1D(hist_name, f"Event {i_ev};Drift Time [ns];Counts", 500, 0, 0.5)

            for p_obj in prop_objs:
                if p_obj.getType() != 255: continue
                mcp_ptr = p_obj.getMCParticle()
                if not mcp_ptr or mcp_ptr.GetUniqueID() != target_mcp.particle_id: continue
                dt_ns = p_obj.getGlobalTime() * 1e-3
                dts.append(p_obj.getGlobalTime())
                self.histograms["drift_time_spectrum"].Fill(dt_ns)
                if self.save_event_hists: self.event_histograms[hist_name].Fill(dt_ns)
                if self.fill_3d:
                    p_pos = p_obj.getLocalPosition()
                    buffer["drift_time_map_xyz"].append((p_pos.X()*1000 - pix_center[0], p_pos.Y()*1000 - pix_center[1], p_pos.Z()*1000 - pix_center[2], dt_ns))

            if dts:
                dts.sort()
                for f in self.drift_fractions:
                    idx = int(len(dts) * (f / 100.0))
                    if idx >= len(dts): idx = len(dts) - 1
                    tf = dts[idx] * 1e-3
                    buffer[f"electron_driftTime_{f}p"].append(tf)
                    buffer[f"inPixel_electron_driftTime_{f}p"].append((in_x, in_y, tf))
                    if f == 90:
                        buffer["cluster_charge_vs_drift_time_90p"].append((q_sum, tf))
                        buffer["cluster_size_vs_drift_time_90p"].append((cl_size, tf))

            if (i_ev + 1) % batch_size == 0: self._fill_histograms_from_buffer(buffer)
            if (i_ev + 1) % events_per_chunk == 0: gc.collect()

        self._fill_histograms_from_buffer(buffer)
        gc.collect()

    def finalize(self):
        for i, name in enumerate(self.counter_names):
            self.histograms["counters"].SetBinContent(i + 1, self.counters.get(name, 0))

        output_file = ROOT.TFile(self.config.get("output_file_name"), "RECREATE")
        dirs = {
            "inPixel": output_file.mkdir("inPixel"),
            "DriftTime": output_file.mkdir("DriftTime"),
        }
        if self.save_per_size: dirs["PerSize"] = output_file.mkdir("PerClusterSize") #
        if self.save_event_hists: dirs["EventDist"] = output_file.mkdir("EventDriftDistributions")
        
        for name, histo in self.histograms.items():
            if name.startswith("inPixel_"): dirs["inPixel"].cd()
            elif any(x in name for x in ["drift_time", "driftTime"]): dirs["DriftTime"].cd()
            elif "_size_" in name:
                if self.save_per_size: dirs["PerSize"].cd() #
                else: continue
            else: output_file.cd()
            histo.Write()
            
        if self.save_event_hists:
            dirs["EventDist"].cd()
            for h_event in self.event_histograms.values(): h_event.Write()
        output_file.Close()

def run_worker(args):
    start_index, end_index, worker_id, config = args
    config["output_file_name"] = f"analysis_py_part_{worker_id}.root"
    config["start_entry"] = start_index
    config["end_entry"] = end_index
    config["worker_id"] = worker_id
    try:
        analyzer = AnalysisPixelModule(config)
        analyzer.run_analysis(events_per_chunk=2000)
        analyzer.finalize()
        return config["output_file_name"]
    except Exception:
        traceback.print_exc()
        return None

if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True)
    parser.add_argument("-o", "--output", default="analysis_py.root", type=str)
    parser.add_argument("-st", "--seed_threshold", type=int, default=1000)
    parser.add_argument("-nt", "--neighbor_threshold", type=int, default=500)
    parser.add_argument("-p", "--pixel_pitch", type=float, default=22.5)
    parser.add_argument("-n", "--name", type=str, default="CE65")
    parser.add_argument("--one_bit", action="store_true")
    parser.add_argument("--fill3D", action="store_true")
    parser.add_argument("--save_event_hists", action="store_true")
    parser.add_argument("--save_per_size", action="store_true", help="Save histograms split by cluster size") #
    parser.add_argument("-j", "--cores", type=int, default=8)

    args = parser.parse_args()
    config = {
        "file_name": args.input_file,
        "output_file_name": args.output,
        "detector_name": args.name,
        "pixel_pitch_x": args.pixel_pitch,
        "pixel_pitch_y": args.pixel_pitch,
        "seed_threshold": args.seed_threshold,
        "neighbor_threshold": args.neighbor_threshold,
        "max_cluster_charge": 20000,
        "one_bit": args.one_bit,
        "fill_3d": args.fill3D,
        "save_event_hists": args.save_event_hists,
        "save_per_size": args.save_per_size, #
    }

    temp_file = ROOT.TFile.Open(config["file_name"])
    total_entries = temp_file.Get("PixelHit").GetEntries()
    temp_file.Close()

    num_cores = min(args.cores, total_entries)
    chunk_size = math.ceil(total_entries / num_cores)
    tasks = [(i*chunk_size, min((i+1)*chunk_size, total_entries), i, config.copy()) for i in range(num_cores)]

    with multiprocessing.Pool(processes=num_cores) as pool:
        partial_files = [res for res in pool.map(run_worker, tasks) if res]

    if partial_files:
        os.system(f"hadd -f {args.output} {' '.join(partial_files)}")
        for f in partial_files: os.remove(f)

    print(f"Finished in {time.time() - start_time:.22f}s.")