import ROOT
import logging
import os
import glob
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import json
import math

# --- Langau関数の定義 (C++マクロの埋め込み) ---
langau_code = """
#include <TMath.h>
#include <TF1.h>

Double_t langaufun(Double_t *x, Double_t *par) {
   Double_t invsq2pi = 0.3989422804014;
   Double_t mpshift  = -0.22278298;
   Double_t np = 100.0;
   Double_t sc =   5.0;
   Double_t xx;
   Double_t mpc;
   Double_t fland;
   Double_t sum = 0.0;
   Double_t xlow,xupp;
   Double_t step;
   Double_t i;

   mpc = par[1] - mpshift * par[0];
   xlow = x[0] - sc * par[3];
   xupp = x[0] + sc * par[3];
   step = (xupp-xlow) / np;

   for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
   }
   return (par[2] * step * sum * invsq2pi / par[3]);
}
"""
ROOT.gInterpreter.Declare(langau_code)

# --- ロギング設定 ---
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s | %(levelname)-7s | %(name)-15s | %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger("plot_Check")

# --- ROOTグローバル設定 ---
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# --- ROOTカラーのインポート ---
from ROOT import kBlue, kAzure, kRed, kPink, kGreen, kOrange, kBlack, kCyan, kGray, kMagenta

# --- DataConfig ---
@dataclass
class DataConfig:
    name: str
    fileKey: str
    scale: float
    color: int
    colorAlt: int
    markerStyle: int = 20
    markerStyleAlt: int = 24
    isMerged: bool = False
    doRebin: bool = True 
    rebinValues: Dict[str, int] = field(default_factory=lambda: {"default": 1})
    histNames: Dict[str, str] = field(default_factory=dict)
    rebinType: str = "default"

def parse_root_color(color_str: str) -> int:
    if isinstance(color_str, int): return color_str
    color_str = color_str.replace(" ", "")
    try:
        if "+" in color_str:
            parts = color_str.split("+")
            return getattr(ROOT, parts[0]) + int(parts[1])
        elif "-" in color_str:
            parts = color_str.split("-")
            return getattr(ROOT, parts[0]) - int(parts[1])
        else:
            return getattr(ROOT, color_str)
    except:
        return ROOT.kBlack

def load_comparisons_config(json_path: str) -> List[Tuple[List[str], str]]:
    """Loads plot comparison combinations from a JSON file."""
    if not os.path.exists(json_path):
        log.error(f"Comparison config file not found: {json_path}")
        return []
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Converts list of dicts to list of tuples (keys, suffix)
    return [(item["keys"], item["suffix"]) for item in data]

def load_data_config(json_path: str) -> Dict[str, DataConfig]:
    if not os.path.exists(json_path):
        log.error(f"Config file not found: {json_path}")
        return {}
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    configs = {}
    for key, val in data.items():
        if "color" in val: val["color"] = parse_root_color(val["color"])
        if "colorAlt" in val: val["colorAlt"] = parse_root_color(val["colorAlt"])
        try:
            configs[key] = DataConfig(**val)
        except TypeError as e:
            log.error(f"Error creating DataConfig for '{key}': {e}")
    return configs

def get_merged_object(base_file_path, hist_name):
    file_list = glob.glob(f"{base_file_path}_*.root")
    if not file_list: return None
    merged_hist = None
    for file_path in file_list:
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie(): continue
        h_obj = f.Get(hist_name)
        if not h_obj or not isinstance(h_obj, ROOT.TH1):
             f.Close(); continue
        if merged_hist is None:
            merged_hist = h_obj.Clone(f"{h_obj.GetName()}_merged")
            merged_hist.SetDirectory(0)
        else:
            merged_hist.Add(h_obj)
        f.Close()
    return merged_hist

class plot_Check:
    def __init__(self, config_file="simulation_data.json"):
        self.m_fileCache = {} 
        self.config_file = config_file
        log.info("plot_Check::__init__ | plot_Check object is created.")
        
        self.output_dirs = {
            "1D": "./plot/1D_Comparison",
            "2D": "./plot/2D_Maps",
            "Path": "./plot/PathAnalysis",
            "SeedCS": "./plot/1D_Comparison/seed_charge_cs"
        }
        for d in self.output_dirs.values():
            os.makedirs(d, exist_ok=True)

    def __del__(self):
        self.cleanupResources()

    def cleanupResources(self):
        for file_name, file in self.m_fileCache.items():
            if file: file.Close()
        self.m_fileCache.clear()

    def openFile(self, fileName):
        if fileName not in self.m_fileCache:
            file = ROOT.TFile.Open(fileName)
            if not file or file.IsZombie():
                self.m_fileCache[fileName] = None
            else:
                self.m_fileCache[fileName] = file
        return self.m_fileCache.get(fileName)

    def setHistStyle(self, hist, color, marker, alpha=0.15):
        if not hist: return
        hist.SetMarkerColor(color)
        hist.SetLineColor(color)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(0.8)
        hist.SetFillColorAlpha(color, alpha)
        hist.SetLineWidth(2)

    def setAxisStyle(self, hist):
        if not hist: return
        label_size = 0.04
        title_size = 0.05
        xaxis = hist.GetXaxis()
        xaxis.SetLabelSize(label_size); xaxis.SetTitleSize(title_size); xaxis.SetTitleOffset(0.9)
        yaxis = hist.GetYaxis()
        yaxis.SetLabelSize(label_size); yaxis.SetTitleSize(title_size); yaxis.SetTitleOffset(1)

    def scalingHistogram(self, hist, scale_factor, title):
        if not hist: return None
        if scale_factor == 1.0:
            hist.SetTitle(title)
            return hist
        xAxis = hist.GetXaxis()
        xAxis.Set(xAxis.GetNbins(), xAxis.GetXmin() * scale_factor, xAxis.GetXmax() * scale_factor)
        hist.SetTitle(title)
        for i in range(0, xAxis.GetNbins() + 2):
            if hist.GetBinContent(i) < 0: hist.SetBinContent(i, 0)
        return hist

    def getScaledHist(self, fileKey, histName, scaleFactor, title, isMerged=False):
        hist = None
        if isMerged:
            hist = get_merged_object(fileKey, histName)
        else:
            f = self.openFile(fileKey)
            if f:
                h_obj = f.Get(histName)
                if h_obj:
                    hist = h_obj.Clone(f"{h_obj.GetName()}_clone")
                    hist.SetDirectory(0)
        if not hist: return None

        is_neighbor_hist = ("neighbor" in histName.lower())
        is_exp_data = ("sps" in fileKey)
        # Check if it is strictly a 1D histogram (TH1 but not TH2)
        is_2d_hist = hist.InheritsFrom("TH2")

        # --- Change: Added check 'not is_2d_hist' to prevent error ---
        if is_neighbor_hist and is_exp_data and not is_2d_hist:
            threshold_adc = 100.0
            for i in range(0, hist.GetNbinsX() + 2):
                if hist.GetBinCenter(i) < threshold_adc:
                    hist.SetBinContent(i, 0)
                    hist.SetBinError(i, 0)
            log.info(f"getScaledHist | Applied cut < {threshold_adc} (ADC) for {histName}.")

        # "sim" がファイルキーに含まれ、かつ2Dヒストグラムの場合
        if "sim" in fileKey.lower() and hist.InheritsFrom("TH2"):
            # ビン数が多すぎる（例: 100以上）場合はリビニングする
            # Simの解像度が高すぎるため、実験データに合わせて粗くする(例: 4x4結合)
            if hist.GetNbinsX() > 100: 
                hist.Rebin2D(4, 4)
                log.info(f"getScaledHist | Rebinning simulation 2D hist {histName} by 4x4 to smooth noise.")
            elif hist.GetNbinsX() > 50:
                hist.Rebin2D(2, 2)

        return self.scalingHistogram(hist, scaleFactor, title)

    def fitLandau(self, hist, xMin, xMax, color):
        func_name = f"fit_{hist.GetName()}_{ROOT.gRandom.Integer(10000)}"
        mean = hist.GetMean()
        sigma_raw = hist.GetRMS()
        max_bin = hist.GetMaximumBin()
        peak_x = hist.GetXaxis().GetBinCenter(max_bin)
        
        fit_min_limit = peak_x - (sigma_raw * 0.8) 
        fit_max_limit = peak_x + (sigma_raw * 3.0)
        if fit_min_limit < xMin: fit_min_limit = xMin
        if fit_max_limit > xMax: fit_max_limit = xMax
        if fit_min_limit >= peak_x: fit_min_limit = peak_x - 100

        func = ROOT.TF1(func_name, "landau", xMin, xMax)
        func.SetParameters(hist.GetMaximum() * 5.0, peak_x, sigma_raw * 0.15)
        func.SetParLimits(1, peak_x - sigma_raw, peak_x + sigma_raw)

        func.SetRange(fit_min_limit, fit_max_limit)
        hist.Fit(func, "Q R 0 N")
        func.SetRange(xMin, xMax) 
        func.SetLineColor(color); func.SetLineStyle(2); func.SetLineWidth(2)
        return func
    
    def fitLangau(self, hist, xMin, xMax, color):
        func_name = f"fit_langau_{hist.GetName()}_{ROOT.gRandom.Integer(10000)}"
        max_bin = hist.GetMaximumBin()
        peak_x = hist.GetXaxis().GetBinCenter(max_bin)
        area = hist.Integral() * hist.GetBinWidth(1)
        rms = hist.GetRMS()

        func = ROOT.TF1(func_name, ROOT.langaufun, xMin, xMax, 4)
        func.SetParNames("LandauWidth", "MPV", "Area", "GausSigma")
        func.SetParameters(rms * 0.15, peak_x, area, rms * 0.10)
        func.SetParLimits(0, 0.1, rms)
        func.SetParLimits(1, peak_x - rms, peak_x + rms)
        func.SetParLimits(2, area * 0.1, area * 10.0)
        func.SetParLimits(3, 1.0, rms)
        func.SetNpx(500)

        fit_range_min = peak_x - rms
        fit_range_max = peak_x + 3 * rms
        if fit_range_min < xMin: fit_range_min = xMin
        if fit_range_max > xMax: fit_range_max = xMax
        
        func.SetRange(fit_range_min, fit_range_max)
        hist.Fit(func, "Q R 0 N")
        func.SetRange(xMin, xMax)
        func.SetLineColor(color); func.SetLineStyle(2); func.SetLineWidth(2)
        hist.Fit(func, "Q R B M 0 N") 
        return func

    # def fitLangau(self, hist, xMin, xMax, color):
    #     """
    #     Final refined Langau fit. 
    #     Focuses on a very narrow asymmetric window to prioritize the peak position 
    #     over the tail shape, especially for simulation data.
    #     """
    #     func_name = f"fit_langau_{hist.GetName()}_{ROOT.gRandom.Integer(10000)}"
    #     max_bin = hist.GetMaximumBin()
    #     peak_x = hist.GetXaxis().GetBinCenter(max_bin)
    #     area = hist.Integral() * hist.GetBinWidth(1)

    #     func = ROOT.TF1(func_name, ROOT.langaufun, xMin, xMax, 4)
    #     func.SetParNames("LandauWidth", "MPV", "Area", "GausSigma")
        
    #     # --- Strict Parameter Initialization ---
    #     # Narrow Landau and Gaussian widths to force a sharper peak
    #     func.SetParameters(30.0, peak_x, area * 0.4, 25.0)
        
    #     # --- Asymmetric Range Adjustment ---
    #     # Simulation has a larger tail, so we cut the right side earlier.
    #     # Experiment (~600e) and Simulation might need slightly different windows, 
    #     # but 180e offset on the right is generally safe to avoid the shoulder.
    #     fit_range_min = peak_x - 120.0
    #     fit_range_max = peak_x + 180.0 

    #     if fit_range_min < xMin: fit_range_min = xMin
    #     if fit_range_max > xMax: fit_range_max = xMax
        
    #     # --- Strict Limits ---
    #     func.SetParLimits(0, 5.0, 150.0)             # Landau Width
    #     func.SetParLimits(1, peak_x - 50, peak_x + 50) # Very tight MPV constraint
    #     func.SetParLimits(3, 5.0, 100.0)             # Gaus Sigma
        
    #     func.SetNpx(500)

    #     # Fit with "L" (Likelihood) and "I" (Integral) for more robust peak finding
    #     # "M" improves the fit result by searching for better parameter values
    #     hist.Fit(func, "Q R L I 0 N") 
        
    #     func.SetLineColor(color)
    #     func.SetLineStyle(2)
    #     func.SetLineWidth(2)
        
    #     # Final pass with Minuit optimization
    #     hist.Fit(func, "Q R L I B M 0 N") 
        
    #     func.SetRange(xMin, xMax)
        
    #     return func
    
    def plot2DComparison(self, c, l, datasets, suffix, histKey, histLegendSuffix, outName, plotTitle, doProfile=False):
        c.Clear(); l.Clear()
        ROOT.gStyle.SetPalette(ROOT.kBird)
        valid_datasets = [ds for ds in datasets if ds.histNames.get(histKey)]
        if not valid_datasets: return

        if doProfile:
            self.setCanvasMargins(c)
            profiles = []; max_y = 0.0
            for ds in valid_datasets:
                h2 = self.getScaledHist(ds.fileKey, ds.histNames.get(histKey), 1.0, "", ds.isMerged)
                if not h2 or not h2.InheritsFrom("TH2"): continue
                
                prof = h2.ProfileX(f"{h2.GetName()}_prof_{ds.name}_{suffix}")
                self.setHistStyle(prof, ds.color, ds.markerStyle)
                prof.SetLineWidth(2)
                profiles.append((prof, ds.name))
                
                if prof.GetMaximum() > max_y: max_y = prof.GetMaximum()

            if not profiles: return
            for i, (prof, name) in enumerate(profiles):
                opt = "PE" if i == 0 else "same PE"
                if i == 0:
                    prof.SetTitle(f"{plotTitle} (Profile)")
                    prof.GetYaxis().SetTitle(f"Average {prof.GetYaxis().GetTitle()}")
                    prof.GetYaxis().SetRangeUser(0, max_y * 1.5)
                    self.setAxisStyle(prof)
                prof.Draw(opt)
                l.AddEntry(prof, f"{name}", "pl")
            
            l.Draw(); self.drawTitle(f"{plotTitle} (Trend)")
            c.SaveAs(f"{self.output_dirs['2D']}/{outName}_ProfileCompare_{suffix}.pdf")
        else:
            n_plots = len(valid_datasets)
            n_cols = math.ceil(math.sqrt(n_plots))
            n_rows = math.ceil(n_plots / n_cols)
            c.Divide(n_cols, n_rows)
            kept_hists = [] 
            for i, ds in enumerate(valid_datasets):
                c.cd(i + 1)
                ROOT.gPad.SetLeftMargin(0.12); ROOT.gPad.SetRightMargin(0.12)
                h2 = self.getScaledHist(ds.fileKey, ds.histNames.get(histKey), 1.0, "", ds.isMerged)
                if h2 and h2.InheritsFrom("TH2"):
                    h2.SetTitle(f"{ds.name}"); self.setAxisStyle(h2)
                    h2.SetStats(0); h2.Draw("COLZ")
                    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.04)
                    t.DrawLatex(0.15, 0.85, f"{ds.name}")
                    kept_hists.append(h2); kept_hists.append(t)
            c.SaveAs(f"{self.output_dirs['2D']}/{outName}_Maps_{suffix}.pdf")

    def run_plotCheck(self, comp_config_file="comparisons_config.json"):
        log.info("run_plotCheck | Start")
        data_entries = load_data_config(self.config_file)
        if not data_entries: return

        charge_rebin_exp_defaults = {"clusterCharge": 20, "seedCharge": 20, "neighborChargeSum": 10, "seedChargeBySize": 20, "clusterSize": 1, "driftTime": 1, "default": 1}
        charge_rebin_sim_defaults = {"clusterCharge": 4, "seedCharge": 4, "neighborChargeSum": 1, "seedChargeBySize": 4, "clusterSize": 1, "driftTime": 1, "electron_driftTime_90p": 100, "default": 1}
        for key, config in data_entries.items():
            if "sim" in key.lower() or "sim" in config.name.lower(): defaults = charge_rebin_sim_defaults.copy()
            else: defaults = charge_rebin_exp_defaults.copy()
            defaults.update(config.rebinValues)
            config.rebinValues = defaults

        # comparisons = load_comparisons_config(comp_config_file)
        # if not comparisons:
        #     log.warning("No comparison sets found. Check your JSON path.")
        #     return

        all_comparisons = load_comparisons_config(comp_config_file)
        if not all_comparisons:
            log.warning(f"No comparison sets found in {comp_config_file}.")
            return

        #target_keywords = ["doping", "efield"] 
        target_keywords = ["p15_STD_trap", "ClusterSize", "ClusterSeedCharge", "ClusterCharge"] 
        comparisons = [
            (keys, suffix) for keys, suffix in all_comparisons 
            if any(key in suffix for key in target_keywords)
        ]

        if not comparisons:
            log.info("No matching suffixes found for the given keywords.")
            return

        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        legend = ROOT.TLegend(0.45, 0.55, 0.9, 0.8)
        legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.03)

        for keys, suffix in comparisons:
            datasets = [data_entries[k] for k in keys if k in data_entries]
            if not datasets: continue
            log.info(f"Processing set: {keys}")
            plot_title = " vs ".join([d.name for d in datasets])
            if len(plot_title) > 1: plot_title = f"{suffix}"

            self.plotAllCharge(canvas, legend, datasets, suffix, plot_title)
            self.plotClusterSize(canvas, legend, datasets, suffix, plot_title)
            self.plotComparison(canvas, legend, datasets, suffix, "seedCharge", "seed", "clusterSeedCharge_comp", plot_title, 0, 4000, False) 
            self.plotComparison(canvas, legend, datasets, suffix, "clusterCharge", "cluster", "clusterCharge_comp", plot_title, 0, 3000, True)
            self.plotComparison(canvas, legend, datasets, suffix, "neighborChargeSum", "neighbor sum", "clusterNeighborChargeSum_comp", plot_title, 20, 2000, False)
            self.plotComparison(canvas, legend, datasets, suffix, "driftTime", "drift", "driftTime_comp", plot_title, 0, 5, False)
            self.plotComparison(canvas, legend, datasets, suffix, "electron_driftTime_90p", "drift 90%", "electronDriftTime90p_comp", plot_title, 0, 102, False)
            self.plotSeedChargeByCS(canvas, legend, datasets, suffix, plot_title)

            map_plots_config = [
                # Spatial Maps
                ("map_seedCharge", "Seed Charge Map", "2D_Map_SeedCharge", False, True, "Average Seed Charge [e]"),
                ("map_clusterSize", "Cluster Size Map", "2D_Map_ClusterSize", False, True, "Average Cluster Size"),
                ("map_neighborChargeSum", "Neighbor Sum Map", "2D_Map_NeighborSum", False, True, "Average Neighbor Charge [e]"),
                
                # Correlation Maps (Y axis titles for these are handled by profile logic, but added for consistency)
                ("map_clusterVSCSize", "Cluster Charge vs Size", "2D_Corr_ClusterVSSize", True, False, "Cluster Charge [e]"),
                ("map_seedVSCSize", "Seed Charge vs Size", "2D_Corr_SeedVSSize", True, False, "Seed Charge [e]"),
                ("map_seedVSneighborChargeSum", "Seed vs Neighbor Sum", "2D_Corr_SeedVSNeighbor", False, False, "Neighbor Sum [e]"),
            ]

            for hist_key, leg_suffix, out_name, do_profile, do_path_analysis, y_title in map_plots_config:
                if any(ds.histNames.get(hist_key) for ds in datasets):
                    if do_profile:
                        self.plot2DComparison(canvas, legend, datasets, suffix, hist_key, leg_suffix, out_name, plot_title, True)
                    elif do_path_analysis:
                        self.plotPathAnalysisCombined(canvas, legend, datasets, hist_key, suffix, plot_title, yAxisTitle=y_title, outName=out_name)

                    c_div = ROOT.TCanvas("c_div", "c_div", 1200, 800)
                    self.plot2DComparison(c_div, legend, datasets, suffix, hist_key, leg_suffix, out_name, plot_title, False)
                    c_div.Close()
                    
                    if do_path_analysis:
                        for ds in datasets:
                            if ds.histNames.get(hist_key):
                                self.plotPathAnalysis(ds, hist_key, suffix)

            canvas.cd()
        
        self.cleanupResources()
        log.info("run_plotCheck | Finished")

    def setCanvasMargins(self, canvas, left=0.14, right=0.06, top=0.10, bottom=0.14):
        if not canvas: return
        canvas.SetLeftMargin(left); canvas.SetRightMargin(right)
        canvas.SetTopMargin(top); canvas.SetBottomMargin(bottom)

    def plotAllCharge(self, c, l, datasets, suffix, plotTitle):
        c.Clear(); l.Clear(); self.setCanvasMargins(c)
        hists_cl = []; hists_sd = []; funcs = []; max_y = 0.0
        for ds in datasets:
            h_cl = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterCharge"), ds.scale, "", ds.isMerged)
            h_sd = self.getScaledHist(ds.fileKey, ds.histNames.get("seedCharge"), ds.scale, "", ds.isMerged)
            if h_cl and h_sd:
                self.setHistStyle(h_cl, ds.color, ds.markerStyle)
                self.setHistStyle(h_sd, ds.colorAlt, ds.markerStyleAlt)
                if ds.doRebin:
                    h_cl.Rebin(ds.rebinValues.get("clusterCharge", 1))
                    h_sd.Rebin(ds.rebinValues.get("seedCharge", 1))
                if h_cl.GetMaximum() > 0: h_cl.Scale(1.0 / h_cl.GetMaximum())
                if h_sd.GetMaximum() > 0: h_sd.Scale(1.0 / h_sd.GetMaximum())
                max_y = max(max_y, h_cl.GetMaximum(), h_sd.GetMaximum())
                funcs.append(self.fitLandau(h_cl, 0, 4000, ds.color))
                hists_cl.append((h_cl, ds.name)); hists_sd.append((h_sd, ds.name))

        if not hists_cl: return
        for i, ((h_cl, name), (h_sd, _)) in enumerate(zip(hists_cl, hists_sd)):
            opt = "PE" if i == 0 else "same PE"
            if i == 0:
                h_cl.SetTitle(";charge [e];counts"); h_cl.GetXaxis().SetRangeUser(0, 4000); h_cl.GetYaxis().SetRangeUser(0, max_y * 1.2)
            h_cl.Draw(opt); h_sd.Draw("same PE"); funcs[i].Draw("same")
            l.AddEntry(h_sd, f"{name}, seed", "pe"); l.AddEntry(h_cl, f"{name}, cluster", "pe")
        l.Draw(); self.drawTitle(plotTitle)
        c.SaveAs(f"{self.output_dirs['1D']}/AllclusterCharge_{suffix}.pdf")

    def plotClusterSize(self, c, l, datasets, suffix, plotTitle):
        c.Clear(); l.Clear(); self.setCanvasMargins(c)
        hists = []; max_y = 0.0
        for ds in datasets:
            h = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterSize"), 1.0, "", ds.isMerged)
            if h:
                self.setHistStyle(h, ds.color, ds.markerStyle, 0.2)
                if h.Integral() > 0: h.Scale(1.0 / h.Integral())
                max_y = max(max_y, h.GetMaximum())
                hists.append((h, ds.name))
        if not hists: return
        for i, (h, name) in enumerate(hists):
            opt_hist = "HIST" if i == 0 else "same HIST"
            if i == 0:
                h.SetTitle(";cluster size;counts"); h.GetXaxis().SetRangeUser(1, 10); h.GetYaxis().SetRangeUser(0, max(max_y * 1.2, 0.1))
            self.setAxisStyle(h); h.Draw(opt_hist); h.Draw("same PE")
            l.AddEntry(h, f"{name} (mean={h.GetMean():.2f})", "pef")
        l.Draw(); self.drawTitle(plotTitle)
        c.SaveAs(f"{self.output_dirs['1D']}/clusterSize_{suffix}.pdf")

    # def plotComparison(self, c, l, datasets, suffix, histKey, histLegendSuffix, outName, plotTitle, xMin, xMax, doLandauFit=False):
    #     c.Clear(); l.Clear(); self.setCanvasMargins(c)
    #     hists = []; funcs = []; max_y = 0.0
    #     title = ";charge [e];counts" if histKey != "driftTime" else ";drift time [ps];counts"
    #     if histKey == "electron_driftTime_90p":
    #         title = ";90% drift time [ps];counts"

    #     for ds in datasets:
    #         current_scale = ds.scale
    #         if histKey in ["driftTime", "electron_driftTime_90p"]: 
    #             current_scale = 1000.0

    #         h = self.getScaledHist(ds.fileKey, ds.histNames.get(histKey), current_scale, "", ds.isMerged)
    #         if h:
    #             self.setHistStyle(h, ds.color, ds.markerStyle)
    #             if ds.doRebin: h.Rebin(ds.rebinValues.get(histKey, 1))
    #             if h.GetMaximum() > 0:
    #                 max_val_limited = 0
    #                 for b in range(1, h.GetNbinsX()-100):
    #                     content = h.GetBinContent(b)
    #                     if content > max_val_limited:
    #                         max_val_limited = content

    #                 scale_ref = max_val_limited if max_val_limited > 0 else h.GetMaximum()
    #                 h.Scale(1.0 / scale_ref)
    #                 #h.Scale(1.0 / h.GetMaximum())
    #             max_y = max(max_y, h.GetMaximum())
    #             funcs.append(self.fitLangau(h, xMin, xMax, ds.color) if doLandauFit else None)
    #             hists.append((h, ds.name))

    #     if not hists: return
    #     legend_items = []
    #     for i, (h, name) in enumerate(hists):
    #         opt = "E3" if i == 0 else "E3 SAME"
    #         if i == 0:
    #             h.SetTitle(title); h.GetXaxis().SetRangeUser(xMin, xMax); h.GetYaxis().SetRangeUser(0, max_y * 1.3)
    #             self.setAxisStyle(h)
    #         h.Draw(opt)
    #         legText = f"{name}, {histLegendSuffix}"
    #         if funcs[i]:
    #             funcs[i].Draw("same"); legText += f" (MPV={funcs[i].GetParameter(1):.0f})"
    #         #l.AddEntry(h, legText, "PF")
    #         h_leg = h.Clone(f"{h.GetName()}_leg")
    #         h_leg.SetLineWidth(0)
    #         legend_items.append(h_leg)
    #         l.AddEntry(h_leg, legText, "pf")
    #     l.Draw(); self.drawTitle(plotTitle)
    #     c.SaveAs(f"{self.output_dirs['1D']}/{outName}_{suffix}.pdf")

    def plotComparison(self, c, l, datasets, suffix, histKey, histLegendSuffix, outName, plotTitle, xMin, xMax, doLandauFit=False):
        c.Clear(); l.Clear(); self.setCanvasMargins(c)
        hists = []; funcs = []; max_y = 0.0
        # To prevent segfault, keep references to cloned histograms
        kept_legend_hists = [] 
        
        title = ";charge [e];counts" if "Charge" in histKey else ";90% drift time [ps];counts"

        for ds in datasets:
            current_scale = ds.scale
            if histKey in ["driftTime", "electron_driftTime_90p"]: 
                current_scale = 1000.0

            target_hist_name = ds.histNames.get(histKey)
            if not target_hist_name: continue 

            h = self.getScaledHist(ds.fileKey, target_hist_name, current_scale, "", ds.isMerged)
            if h:
                self.setHistStyle(h, ds.color, ds.markerStyle)
                if ds.doRebin: h.Rebin(ds.rebinValues.get(histKey, 1))
                
                # --- Start of modified scaling logic ---
                # If it's a drift time plot, find the peak within a manually defined range
                if "driftTime" in histKey:
                    # Adjust these values based on where your physical peak usually resides
                    search_x_min = 0.0
                    search_x_max = 90.0 
                    
                    b_start = h.FindBin(search_x_min)
                    b_end = h.FindBin(search_x_max)
                    
                    peak_val = 0.0
                    for i in range(b_start, b_end + 1):
                        val = h.GetBinContent(i)
                        if val > peak_val:
                            peak_val = val
                    
                    scale_ref = peak_val if peak_val > 0 else h.GetMaximum()
                else:
                    scale_ref = h.GetMaximum()

                if scale_ref > 0:
                    h.Scale(1.0 / scale_ref)
                
                # Zero out the very last bin (accumulation bin) for a cleaner plot
                h.SetBinContent(h.GetNbinsX(), 0)
                h.SetBinError(h.GetNbinsX(), 0)
                # --- End of modified scaling logic ---
                
                max_y = max(max_y, h.GetMaximum())
                funcs.append(self.fitLangau(h, xMin, xMax, ds.color) if doLandauFit else None)
                hists.append((h, ds.name))

        if not hists: return

        for i, (h, name) in enumerate(hists):
            opt = "E3" if i == 0 else "E3 SAME"
            if i == 0:
                h.SetTitle(title)
                h.GetXaxis().SetRangeUser(xMin, xMax)
                h.GetYaxis().SetRangeUser(0, 1.3) 
                self.setAxisStyle(h)
            
            h.Draw(opt)
            h.Draw("same PE")
            
            legText = f"{name}, {histLegendSuffix}"
            if funcs[i]:
                funcs[i].Draw("same")
            
            h_leg = h.Clone(f"{h.GetName()}_leg")
            h_leg.SetDirectory(0)
            h_leg.SetLineWidth(0)
            kept_legend_hists.append(h_leg)
            l.AddEntry(h_leg, legText, "pf")
            
        l.Draw(); self.drawTitle(plotTitle)
        c.SaveAs(f"{self.output_dirs['1D']}/{outName}_{suffix}.pdf")

    def plotSeedChargeByCS(self, c, l, datasets, suffix, plotTitleBase):
        valid_datasets = [ds for ds in datasets if ds.histNames.get("seedChargeBySize")]
        if not valid_datasets: return
        for cs in range(1, 7):
            c.Clear(); l.Clear(); self.setCanvasMargins(c)
            hists = []; max_y = 0.0
            for ds in valid_datasets:
                baseName = ds.histNames.get("seedChargeBySize")
                targetName = f"PerClusterSize/clsize_{cs}/seed_charge_size_{cs}" if "Sim" in ds.name else f"{baseName}{cs}"
                h = self.getScaledHist(ds.fileKey, targetName, ds.scale, "", ds.isMerged)
                if h:
                    self.setHistStyle(h, ds.color, ds.markerStyle)
                    if ds.doRebin: h.Rebin(ds.rebinValues.get("seedChargeBySize", 1))
                    if h.GetMaximum() > 0: h.Scale(1.0 / h.GetMaximum())
                    max_y = max(max_y, h.GetMaximum())
                    hists.append((h, ds.name))
            if not hists: continue
            for i, (h, name) in enumerate(hists):
                opt = "PE" if i == 0 else "same PE"
                if i == 0:
                    h.SetTitle(f"Seed Charge (Size {cs});charge [e];counts"); h.GetXaxis().SetRangeUser(0, 4000); h.GetYaxis().SetRangeUser(0, max_y * 1.2)
                h.Draw(opt); l.AddEntry(h, f"{name}, seed", "pe")
            l.Draw(); self.drawTitle(f"{plotTitleBase} size {cs}")
            c.SaveAs(f"{self.output_dirs['SeedCS']}/clusterSeedCharge_cs{cs}_{suffix}.pdf")

    def drawTitle(self, text):
        t = ROOT.TLatex(); t.SetTextAlign(12); t.SetTextSize(0.05); t.DrawLatexNDC(0.5, 0.85, text)

    def extractDataAlongPath(self, hist, x1, y1, x2, y2, dist_offset, n_steps=None): # n_stepsを引数から外すかデフォルトNoneに
        if not hist: return None
        x_vals, y_vals, x_errs, y_errs = [], [], [], []
        segment_len = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        # --- 変更: ステップ数を距離に基づいて決定 (約 0.8 um 刻み) ---
        step_size = 0.8  # 実験データのプロットの粒度に合わせる
        if n_steps is None:
            n_steps = int(segment_len / step_size)
            if n_steps < 1: n_steps = 1
        # -------------------------------------------------------------
        
        for i in range(n_steps + 1):
            t = i / float(n_steps)
            curr_x = x1 + t * (x2 - x1)
            curr_y = y1 + t * (y2 - y1)
            
            # Interpolateを使用
            val = hist.Interpolate(curr_x, curr_y)
            
            bin_idx = hist.FindBin(curr_x, curr_y)
            err = hist.GetBinError(bin_idx)

            x_vals.append(dist_offset + t * segment_len)
            y_vals.append(val)
            x_errs.append(0.0)
            y_errs.append(err)
            
        if not x_vals: return None
        import array
        return ROOT.TGraphErrors(len(x_vals), array.array('d', x_vals), array.array('d', y_vals), array.array('d', x_errs), array.array('d', y_errs))
    
    
    def plotPathAnalysis(self, ds, histKey, suffix):
        histName = ds.histNames.get(histKey)
        if not histName: return None
        h2 = self.getScaledHist(ds.fileKey, histName, 1.0, "", ds.isMerged)
        if not h2: return None
        
        c_path = ROOT.TCanvas(f"c_path_{ds.name}", "Path Analysis", 1200, 600)
        c_path.Divide(2, 1)
        
        pad1 = c_path.cd(1); pad1.SetRightMargin(0.15); pad1.SetLeftMargin(0.12)
        h2.SetStats(0); h2.SetTitle(f"{ds.name} 2D Map"); self.setAxisStyle(h2)
        h2.Draw("COLZ")
        
        x_max_axis = h2.GetXaxis().GetXmax(); y_max_axis = h2.GetYaxis().GetXmax()
        inset_x = x_max_axis * 0.1; inset_y = y_max_axis * 0.1
        pA = (0.0, 0.0)
        pB = (0.0, y_max_axis - inset_y)
        pC = (x_max_axis - inset_x, y_max_axis - inset_y)
        pD = (x_max_axis - inset_x, 0.0)
        
        arrow_size = 0.02; arrows = []
        arrows.append(ROOT.TArrow(pA[0], pA[1], pB[0], pB[1], arrow_size, ">"))
        arrows.append(ROOT.TArrow(pB[0], pB[1], pC[0], pC[1], arrow_size, ">"))
        arrows.append(ROOT.TArrow(pC[0], pC[1], pA[0], pA[1], arrow_size, ">"))
        arrows.append(ROOT.TArrow(pA[0], pA[1], pD[0], pD[1], arrow_size, ">"))
        
        colors = [ROOT.kOrange+7, ROOT.kPink-3, ROOT.kRed+1, ROOT.kViolet-2]
        for i, arr in enumerate(arrows):
            arr.SetLineColorAlpha(colors[i], 0.8); arr.SetLineWidth(3); arr.Draw()
            
        latex = ROOT.TLatex(); latex.SetTextSize(0.05)
        latex.DrawLatex(pA[0], pA[1], "A"); latex.DrawLatex(pB[0], pB[1], "B")
        latex.DrawLatex(pC[0], pC[1], "C"); latex.DrawLatex(pD[0], pD[1], "D")

        pad2 = c_path.cd(2); pad2.SetLeftMargin(0.15); pad2.SetGridy()
        
        dist_AB = math.sqrt((pB[0]-pA[0])**2 + (pB[1]-pA[1])**2)
        dist_BC = math.sqrt((pC[0]-pB[0])**2 + (pC[1]-pB[1])**2)
        dist_CA = math.sqrt((pA[0]-pC[0])**2 + (pA[1]-pC[1])**2)

        g_AB = self.extractDataAlongPath(h2, pA[0], pA[1], pB[0], pB[1], 0.0)
        g_BC = self.extractDataAlongPath(h2, pB[0], pB[1], pC[0], pC[1], dist_AB)
        g_CA = self.extractDataAlongPath(h2, pC[0], pC[1], pA[0], pA[1], dist_AB + dist_BC)
        g_AD = self.extractDataAlongPath(h2, pA[0], pA[1], pD[0], pD[1], dist_AB + dist_BC + dist_CA)

        mg = ROOT.TMultiGraph()
        mg.SetTitle(f";Distance along path [um];{h2.GetZaxis().GetTitle()}")
        graphs = [g_AB, g_BC, g_CA, g_AD]
        labels = ["A #rightarrow B", "B #rightarrow C", "C #rightarrow A", "A #rightarrow D"]
        legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9); legend.SetTextSize(0.04)
        
        kept_bands = []

        for i, g in enumerate(graphs):
            if g:
                g.SetLineColor(colors[i]); g.SetLineWidth(2); g.SetMarkerColor(colors[i]); g.SetMarkerStyle(20); g.SetMarkerSize(0.5)
                g_band = g.Clone(); g_band.SetFillColorAlpha(colors[i], 0.3); g_band.SetFillStyle(1001)
                kept_bands.append(g_band)
                mg.Add(g_band, "3"); mg.Add(g, "L"); legend.AddEntry(g, labels[i], "l")
        
        mg.Draw("A"); mg.GetYaxis().SetRangeUser(0, h2.GetMaximum() * 1.2)
        legend.Draw()
        
        line = ROOT.TLine(); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+1)
        boundaries = [0.0, dist_AB, dist_AB+dist_BC, dist_AB+dist_BC+dist_CA]
        bound_labels = ["A", "B", "C", "A"]
        for x, lbl in zip(boundaries, bound_labels):
            line.DrawLine(x, mg.GetYaxis().GetXmin(), x, mg.GetYaxis().GetXmax())
            latex.DrawLatex(x, mg.GetYaxis().GetXmax() * 0.9, lbl)

        c_path.SaveAs(f"{self.output_dirs['Path']}/PathAnalysis_{ds.name}_{suffix}.pdf")
        return mg

    def plotPathAnalysisCombined(self, c, l, datasets, histKey, suffix, plotTitle, yAxisTitle="Value", outName="CombinedPath"):
        c.Clear(); l.Clear(); self.setCanvasMargins(c)
        mg_combined = ROOT.TMultiGraph()
        # Set titles later via mg_combined.GetHistogram() after Draw to ensure it works
        
        has_data = False
        kept_objs = [] # Keep references to prevent garbage collection

        # Temporary list to store boundaries from the first dataset
        boundaries = []
        bound_labels = ["A", "B", "C", "A", "D"]

        for ds in datasets:
            histName = ds.histNames.get(histKey)
            if not histName: continue
            h2 = self.getScaledHist(ds.fileKey, histName, 1.0, "", ds.isMerged)
            if not h2: continue
            
            x_max_axis = h2.GetXaxis().GetXmax(); y_max_axis = h2.GetYaxis().GetXmax()
            inset_x = x_max_axis * 0.1; inset_y = y_max_axis * 0.1
            pA, pB, pC, pD = (0,0), (0, y_max_axis-inset_y), (x_max_axis-inset_x, y_max_axis-inset_y), (x_max_axis-inset_x, 0)
            
            dist_AB = math.sqrt((pB[0]-pA[0])**2 + (pB[1]-pA[1])**2)
            dist_BC = math.sqrt((pC[0]-pB[0])**2 + (pC[1]-pB[1])**2)
            dist_CA = math.sqrt((pA[0]-pC[0])**2 + (pA[1]-pC[1])**2)
            
            if not boundaries:
                boundaries = [0.0, dist_AB, dist_AB+dist_BC, dist_AB+dist_BC+dist_CA, dist_AB+dist_BC+dist_CA+dist_AB]

            g_parts = [
                self.extractDataAlongPath(h2, pA[0], pA[1], pB[0], pB[1], 0.0),
                self.extractDataAlongPath(h2, pB[0], pB[1], pC[0], pC[1], dist_AB),
                self.extractDataAlongPath(h2, pC[0], pC[1], pA[0], pA[1], dist_AB + dist_BC),
                self.extractDataAlongPath(h2, pA[0], pA[1], pD[0], pD[1], dist_AB + dist_BC + dist_CA)
            ]
            
            x_buf, y_buf, ex_buf, ey_buf = [], [], [], []
            for g in g_parts:
                if not g: continue
                for k in range(g.GetN()):
                    x_buf.append(g.GetX()[k]); y_buf.append(g.GetY()[k]); ex_buf.append(0); ey_buf.append(g.GetEY()[k])
            
            if not x_buf: continue
            import array
            # g_total = ROOT.TGraphErrors(len(x_buf), array.array('d', x_buf), array.array('d', y_buf), array.array('d', ex_buf), array.array('d', ey_buf))
            # g_total.SetLineColor(ds.color); g_total.SetLineWidth(2); g_total.SetMarkerStyle(ds.markerStyle); g_total.SetMarkerColor(ds.color); g_total.SetMarkerSize(0.6)
            g_total = ROOT.TGraph(len(x_buf), array.array('d', x_buf), array.array('d', y_buf))
            g_total.SetLineColor(ds.color); g_total.SetLineWidth(2)
            g_total.SetMarkerStyle(ds.markerStyle); g_total.SetMarkerColor(ds.color); g_total.SetMarkerSize(0.7)
            
            # Band style (similar to E3)
            # #g_band = g_total.Clone(f"band_{ds.name}_{ROOT.gRandom.Integer(1000)}")
            # g_band.SetFillColorAlpha(ds.color, 0.15); g_band.SetFillStyle(1001); g_band.SetLineWidth(0)
            g_band = ROOT.TGraphErrors(len(x_buf), array.array('d', x_buf), array.array('d', y_buf), array.array('d', ex_buf), array.array('d', ey_buf))
            g_band.SetFillColorAlpha(ds.color, 0.2); g_band.SetFillStyle(1001); g_band.SetLineWidth(0); g_band.SetMarkerSize(0.7); g_band.SetMarkerStyle(ds.markerStyle); g_band.SetMarkerColor(ds.color)
            
            mg_combined.Add(g_band, "3")
            mg_combined.Add(g_total, "PL")
            
            # Legend entry (matching plotComparison style)
            l.AddEntry(g_band, ds.name, "plf")
            
            kept_objs.extend([g_total, g_band])
            has_data = True

        if not has_data: return

        mg_combined.Draw("A")
        mg_combined.SetTitle(f"{plotTitle};Distance along path [#mum];{yAxisTitle}")
        self.setAxisStyle(mg_combined.GetHistogram())
        
        # Add boundary lines and labels
        line = ROOT.TLine(); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+1)
        latex = ROOT.TLatex(); latex.SetTextSize(0.04); latex.SetTextColor(ROOT.kGray+2)
        y_max = mg_combined.GetYaxis().GetXmax()
        for x, lbl in zip(boundaries, bound_labels):
            line.DrawLine(x, mg_combined.GetYaxis().GetXmin(), x, y_max)
            latex.DrawLatex(x + 0.5, y_max * 0.92, lbl)

        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"{self.output_dirs['Path']}/CombinedPathAnalysis_{outName}_{suffix}.pdf")

if __name__ == "__main__":
    plotter = plot_Check("modeling_check_allData.json")
    #plotter = plot_Check("tuned_analysis.json")
    try:
        plotter.run_plotCheck()
    except Exception as e:
        log.critical(f"Error: {e}", exc_info=True)
    finally:
        plotter.cleanupResources()