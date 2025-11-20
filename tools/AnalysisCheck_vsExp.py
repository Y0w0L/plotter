import ROOT
import logging
import os
import glob
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

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
    """単一のデータソース（実験またはシミュレーション）を定義する"""
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


# --- ヘルパー関数 ---
def get_merged_object(base_file_path, hist_name):
    file_list = glob.glob(f"{base_file_path}_*.root")
    if not file_list:
        log.error(f"get_merged_object | No files found for base path: {base_file_path}")
        return None

    merged_hist = None
    for file_path in file_list:
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie(): continue
        
        h_obj = f.Get(hist_name)
        if not h_obj or not isinstance(h_obj, ROOT.TH1):
             f.Close()
             continue

        if merged_hist is None:
            merged_hist = h_obj.Clone(f"{h_obj.GetName()}_merged")
            merged_hist.SetDirectory(0)
        else:
            merged_hist.Add(h_obj)
        f.Close()
    return merged_hist


class plot_Check:
    def __init__(self):
        self.m_fileCache = {} 
        log.info("plot_Check::__init__ | plot_Check object is created.")
        os.makedirs("./plot/seed_charge_cs", exist_ok=True)

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
                log.error(f"openFile | Failed to open: {fileName}")
                self.m_fileCache[fileName] = None
            else:
                self.m_fileCache[fileName] = file
        return self.m_fileCache.get(fileName)

    def setHistStyle(self, hist, color, marker, alpha=0.5):
        if not hist: return
        hist.SetMarkerColor(color)
        hist.SetLineColor(color)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(0.8)
        hist.SetFillColorAlpha(color, alpha)
        hist.SetLineWidth(2)

    def scalingHistogram(self, hist, scale_factor, title):
        if not hist: return None
        if scale_factor == 1.0:
            hist.SetTitle(title)
            return hist

        xAxis = hist.GetXaxis()
        nbins = xAxis.GetNbins()
        xmin = xAxis.GetXmin()
        xmax = xAxis.GetXmax()

        new_xmin = xmin * scale_factor
        new_xmax = xmax * scale_factor

        xAxis.Set(nbins, new_xmin, new_xmax)
        hist.SetTitle(title)

        for i in range(0, nbins + 2):
            if hist.GetBinContent(i) < 0:
                hist.SetBinContent(i, 0)
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
        return self.scalingHistogram(hist, scaleFactor, title)

    # ★★★ 修正版 fitLandau ★★★
    def fitLandau(self, hist, xMin, xMax, color):
        """
        ヒストグラムに対してLandauフィットを行い、関数オブジェクトを返す。
        フィット範囲の下限をデータの立ち上がりに合わせて自動調整する。
        """
        func_name = f"fit_{hist.GetName()}_{ROOT.gRandom.Integer(10000)}"
        
        maxValue = hist.GetMaximum()
        maxBin = hist.GetMaximumBin()
        peakX = hist.GetXaxis().GetBinCenter(maxBin)
        
        # シグマの初期値を少し大きめに推測 (標準偏差の半分くらい)
        sigma_guess = hist.GetStdDev() * 0.5 

        # ★★★ 自動範囲調整ロジック ★★★
        # ピークの 20% の高さになる最初のビンを探し、そこをフィット開始点とする
        # これにより 0付近の空白データにフィットが引っ張られるのを防ぐ
        threshold_ratio = 0.2
        low_edge_bin = hist.FindFirstBinAbove(maxValue * threshold_ratio)
        fitMin = hist.GetBinCenter(low_edge_bin)
        
        # 計算した開始点がピークより右にあっては困るのでチェック
        if fitMin >= peakX: 
            fitMin = peakX - sigma_guess * 2

        # 指定されたxMinより下にはしない (xMin=0なら0以上になる)
        if fitMin < xMin: fitMin = xMin

        # ログ出力で確認用
        log.info(f"fitLandau | Auto-detected fit range for {hist.GetName()}: {fitMin:.1f} - {xMax:.1f}")

        # Landau関数の定義
        func = ROOT.TF1(func_name, "landau", fitMin, xMax)
        
        # パラメータ初期値: [Constant, MPV, Sigma]
        func.SetParameters(maxValue * 10, peakX, sigma_guess)
        
        # スタイル設定
        func.SetLineColor(color)
        func.SetLineStyle(2) # 破線
        func.SetLineWidth(2)
        
        # フィット実行 (R:指定範囲, Q:静か, 0:描画しない, N:ヒストグラムに保存しない)
        hist.Fit(func, "R Q 0 N")
        
        return func

    def run_plotCheck(self):
        """メイン実行関数"""
        log.info("run_plotCheck | Start")

        # --- データセット定義 ---
        charge_rebin_exp_defaults = {
            "clusterCharge": 20,
            "seedCharge": 20,
            "neighborChargeSum": 20,
            "seedChargeBySize": 20, 
            "clusterSize": 1,
            "default": 1
        }
        charge_rebin_sim_defaults = {
            "clusterCharge": 4,
            "seedCharge": 4,
            "neighborChargeSum": 4,
            "seedChargeBySize": 4, 
            "clusterSize": 1,
            "default": 1
        }

        data_entries = {
            "exp_gap": DataConfig(
                name="Exp GAP (nt200)",
                fileKey="/home/towa/alice3/hist/sps_check/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e",
                scale=0.238,
                color=kBlack, colorAlt=kGray+1,
                markerStyle=20, markerStyleAlt=24,
                isMerged=True,
                doRebin=True,
                rebinValues=charge_rebin_exp_defaults.copy(),
                histNames={
                    "clusterCharge": "AnalysisCE65/CE65_6/cluster/clusterCharge",
                    "seedCharge": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge",
                    "neighborChargeSum": "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum",
                    "clusterSize": "AnalysisCE65/CE65_6/cluster/clusterSize",
                    "seedChargeBySize": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size", 
                }
            ),
            "exp_std": DataConfig(
                name="Exp STD (nt200)",
                fileKey="/home/towa/alice3/hist/sps_check/sps202404_15_std_10V_SeedThd1000e_NeighborThd200e",
                scale=0.240,
                color=kBlack, colorAlt=kGray+1,
                markerStyle=20, markerStyleAlt=24, 
                isMerged=True,
                doRebin=True, 
                rebinValues=charge_rebin_exp_defaults.copy(),
                histNames={
                    "clusterCharge": "AnalysisCE65/CE65_6/cluster/clusterCharge",
                    "seedCharge": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge",
                    "neighborChargeSum": "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum",
                    "clusterSize": "AnalysisCE65/CE65_6/cluster/clusterSize",
                    "seedChargeBySize": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size",
                }
            ),
            "exp_gap_nt600": DataConfig(
                name="Exp GAP (nt600)",
                fileKey="/home/towa/alice3/hist/sps202404/sps202404_15_gap_10V_SeedThd1000e_NeighborThd600e",
                scale=0.238,
                color=kBlack, colorAlt=kGray+1,
                markerStyle=22, markerStyleAlt=26,
                isMerged=True,
                doRebin=True,
                rebinValues=charge_rebin_exp_defaults.copy(),
                histNames={
                    "clusterCharge": "AnalysisCE65/CE65_6/cluster/clusterCharge",
                    "seedCharge": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge",
                    "neighborChargeSum": "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum",
                    "clusterSize": "AnalysisCE65/CE65_6/cluster/clusterSize",
                    "seedChargeBySize": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size", 
                }
            ),
            "exp_std_nt600": DataConfig(
                name="Exp STD nt600",
                fileKey="/home/towa/alice3/hist/sps202404/sps202404_15_std_10V_SeedThd1000e_NeighborThd600e",
                scale=0.240,
                color=kBlack, colorAlt=kGray+1,
                markerStyle=22, markerStyleAlt=26,
                isMerged=True,
                doRebin=True, 
                rebinValues=charge_rebin_exp_defaults.copy(),
                histNames={
                    "clusterCharge": "AnalysisCE65/CE65_6/cluster/clusterCharge",
                    "seedCharge": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge",
                    "neighborChargeSum": "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum",
                    "clusterSize": "AnalysisCE65/CE65_6/cluster/clusterSize",
                    "seedChargeBySize": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size",
                }
            ),
            "sim_gap_mandic_et2": DataConfig(
                name="Sim GAP (mandic et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_mandic.root",
                scale=1,
                color=kAzure+7, colorAlt=kAzure,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_mandic_et2": DataConfig(
                name="Sim STD (mandic et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_mandic.root",
                scale=1,
                color=kAzure+7, colorAlt=kAzure,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_ljubljana_et2": DataConfig(
                name="Sim GAP (ljubljana et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_ljubljana.root",
                scale=1,
                color=kGreen+2, colorAlt=kGreen-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(), 
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_ljubljana_et2": DataConfig(
                name="Sim STD (ljubljana et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_ljubljana.root",
                scale=1,
                color=kGreen+2, colorAlt=kGreen-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(), 
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_dortmund_et2": DataConfig(
                name="Sim GAP (dortmund et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_dortmund.root",
                scale=1,
                color=kMagenta+1, colorAlt=kMagenta-6,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_dortmund_et2": DataConfig(
                name="Sim STD (dortmund et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_dortmund.root",
                scale=1,
                color=kMagenta+1, colorAlt=kMagenta-6,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_cmstracker_et2": DataConfig(
                name="Sim GAP (cmstracker et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_cmstracker.root",
                scale=1,
                color=kRed+1, colorAlt=kRed-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_cmstracker_et2": DataConfig(
                name="Sim STD (cmstracker et2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_cmstracker.root",
                scale=1,
                color=kRed+1, colorAlt=kRed-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_masetti_et2": DataConfig(
                name="Sim GAP (w/oRC nt48e)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_withoutRC.root",
                scale=1,
                #color=kBlue+2, colorAlt=kBlue-2,
                color=kRed+1, colorAlt=kRed-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_masetti_et2": DataConfig(
                name="Sim STD (w/oRC nt48e)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt48e_pip_120GeV_masetti_et2_withoutRC.root",
                scale=1,
                #color=kBlue+2, colorAlt=kBlue-2,
                color=kRed+1, colorAlt=kRed-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_masetti_et2_nt120": DataConfig(
                name="Sim GAP (w/oRC nt120e)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt120e_pip_120GeV_masetti_et2_withoutRC.root",
                scale=1,
                color=kGreen+2, colorAlt=kGreen-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_masetti_et2_nt120": DataConfig(
                name="Sim STD (w/oRC nt120e)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt120e_pip_120GeV_masetti_et2_withoutRC.root",
                scale=1,
                color=kGreen+2, colorAlt=kGreen-4,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_masetti_et2_nt144": DataConfig(
                name="Sim GAP (w/oRC nt144e)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_303k_st240e_nt144e_pip_120GeV_masetti_et2_withoutRC.root",
                scale=1,
                color=kBlue+2, colorAlt=kBlue-2,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_masetti_et2_nt144": DataConfig(
                name="Sim STD (w/oRC nt144e)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_303k_st240e_nt144e_pip_120GeV_masetti_et2_withoutRC.root",
                scale=1,
                color=kBlue+2, colorAlt=kBlue-2,
                markerStyle=20, markerStyleAlt=24,
                isMerged=False,
                doRebin=True, 
                rebinValues=charge_rebin_sim_defaults.copy(),
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
        }

        comparisons = [
            # 1. GAPの Exp vs Sim比較
            (["exp_gap", "exp_gap_nt600", "sim_gap_masetti_et2", "sim_gap_masetti_et2_nt144"], "_gap_default"),
            #(["exp_gap", "sim_gap_masetti_et2"], "_gap_default"),
            #(["exp_gap_nt600", "sim_gap_masetti_et2_nt120", "sim_gap_masetti_et2_nt144"], "_gap_default"),
           
            

            # 2. STDの Exp vs Sim比較
            (["exp_std", "exp_std_nt600", "sim_std_masetti_et2", "sim_std_masetti_et2_nt144"], "_std_default"),
            #(["exp_std", "sim_std_masetti_et2"], "_std_default"),
            #(["exp_std_nt600", "sim_std_masetti_et2_nt120", "sim_std_masetti_et2_nt144"], "_std_default"),
        ]

        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        legend = ROOT.TLegend(0.5, 0.55, 0.9, 0.8)
        legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.025)

        for keys, suffix in comparisons:
            datasets = [data_entries[k] for k in keys if k in data_entries]
            if not datasets: continue

            log.info(f"Processing comparison set: {keys}")
            
            plot_title = " vs ".join([d.name for d in datasets])
            if len(plot_title) > 50: plot_title = f"Comparison {suffix}"

            self.plotAllCharge(canvas, legend, datasets, suffix, plot_title)
            self.plotClusterSize(canvas, legend, datasets, suffix, plot_title)

            self.plotComparison(canvas, legend, datasets, suffix,
                histKey="seedCharge", histLegendSuffix="seed",
                outName="clusterSeedCharge_comp", plotTitle=plot_title,
                xMin=0, xMax=4000,
                doLandauFit=False) 
            
            self.plotComparison(canvas, legend, datasets, suffix,
                histKey="clusterCharge", histLegendSuffix="cluster",
                outName="clusterCharge_comp", plotTitle=plot_title,
                xMin=0, xMax=4000,
                doLandauFit=True) # ClusterChargeはFitする

            self.plotComparison(canvas, legend, datasets, suffix,
                histKey="neighborChargeSum", histLegendSuffix="neighbor sum",
                outName="clusterNeighborChargeSum_comp", plotTitle=plot_title,
                xMin=0, xMax=2000,
                doLandauFit=False)

            self.plotSeedChargeByCS(canvas, legend, datasets, suffix, plot_title)
        
        self.cleanupResources()
        log.info("run_plotCheck | Finished")


    # --- プロット関数群 (変更なし部分は省略せず記載) ---

    def plotAllCharge(self, c, l, datasets: List[DataConfig], suffix: str, plotTitle: str):
        c.Clear(); l.Clear()
        hists_cl = []
        hists_sd = []
        funcs = []
        max_y = 0.0

        for ds in datasets:
            h_cl = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterCharge"), ds.scale, "", ds.isMerged)
            h_sd = self.getScaledHist(ds.fileKey, ds.histNames.get("seedCharge"), ds.scale, "", ds.isMerged)
            
            if h_cl and h_sd:
                self.setHistStyle(h_cl, ds.color, ds.markerStyle)
                self.setHistStyle(h_sd, ds.colorAlt, ds.markerStyleAlt)
                
                if ds.doRebin:
                    rebin_cl = ds.rebinValues.get("clusterCharge", 1)
                    rebin_sd = ds.rebinValues.get("seedCharge", 1)
                    if rebin_cl > 1: h_cl.Rebin(rebin_cl)
                    if rebin_sd > 1: h_sd.Rebin(rebin_sd)
                
                if h_cl.GetMaximum() > 0: h_cl.Scale(1.0 / h_cl.GetMaximum())
                if h_sd.GetMaximum() > 0: h_sd.Scale(1.0 / h_sd.GetMaximum())
                max_y = max(max_y, h_cl.GetMaximum(), h_sd.GetMaximum())
                
                # AllChargeでもフィット線を出したい場合はここでも呼ぶ
                f = self.fitLandau(h_cl, 0, 4000, ds.color)
                funcs.append(f)

                hists_cl.append((h_cl, ds.name))
                hists_sd.append((h_sd, ds.name))

        if not hists_cl: return

        for i, ((h_cl, name), (h_sd, _)) in enumerate(zip(hists_cl, hists_sd)):
            opt = "PE" if i == 0 else "same PE"
            if i == 0:
                h_cl.SetTitle(";charge [e];counts")
                h_cl.GetXaxis().SetRangeUser(0, 4000)
                h_cl.GetYaxis().SetRangeUser(0, max_y * 1.2)
            h_cl.Draw(opt)
            h_sd.Draw("same PE")
            if i < len(funcs): funcs[i].Draw("same")
            
            l.AddEntry(h_sd, f"{name}, seed", "pe")
            l.AddEntry(h_cl, f"{name}, cluster", "pe")

        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"./plot/AllclusterCharge_{suffix}.pdf")

    def plotClusterSize(self, c, l, datasets: List[DataConfig], suffix: str, plotTitle: str):
        c.Clear(); l.Clear()
        hists = []
        max_y = 0.0
        for ds in datasets:
            h = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterSize"), 1.0, "", ds.isMerged) 
            if h:
                self.setHistStyle(h, ds.color, ds.markerStyle, 0.2)
                if h.GetEntries() > 0: h.Scale(1.0 / h.GetEntries())
                max_y = max(max_y, h.GetMaximum())
                hists.append((h, ds.name))

        if not hists: return
        for i, (h, name) in enumerate(hists):
            opt_hist = "HIST" if i == 0 else "same HIST"
            opt_pe = "same PE"
            if i == 0:
                h.SetTitle(";cluster size;counts")
                h.GetXaxis().SetRangeUser(1, 10)
                h.GetYaxis().SetRangeUser(0, max(max_y * 1.2, 0.1)) 
            h.Draw(opt_hist)
            h.Draw(opt_pe) 
            l.AddEntry(h, f"{name} (mean={h.GetMean():.2f})", "pef")
        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"./plot/clusterSize_{suffix}.pdf")

    def plotComparison(self, c, l, datasets: List[DataConfig], suffix: str,
                       histKey: str, histLegendSuffix: str,
                       outName: str, plotTitle: str,
                       xMin: float, xMax: float,
                       doLandauFit: bool = False):
        c.Clear(); l.Clear()
        hists = []
        funcs = [] 
        max_y = 0.0

        for ds in datasets:
            histName = ds.histNames.get(histKey)
            if not histName: continue

            h = self.getScaledHist(ds.fileKey, histName, ds.scale, "", ds.isMerged)
            if h:
                self.setHistStyle(h, ds.color, ds.markerStyle)
                rebin_val = ds.rebinValues.get(histKey, ds.rebinValues.get("default", 1))
                if rebin_val > 1 and ds.doRebin:
                    h.Rebin(rebin_val)
                if h.GetMaximum() > 0: h.Scale(1.0 / h.GetMaximum())
                max_y = max(max_y, h.GetMaximum())

                f = None
                if doLandauFit:
                    # ここでフィットを実行
                    f = self.fitLandau(h, xMin, xMax, ds.color)
                    funcs.append(f)
                else:
                    funcs.append(None)
                hists.append((h, ds.name))

        if not hists: return

        for i, (h, name) in enumerate(hists):
            opt = "PE" if i == 0 else "same PE"
            if i == 0:
                h.SetTitle(";charge [e];counts")
                h.GetXaxis().SetRangeUser(xMin, xMax)
                h.GetYaxis().SetRangeUser(0, max_y * 1.2)
            h.Draw(opt)
            
            legendText = f"{name}, {histLegendSuffix}"
            func = funcs[i]
            if func:
                func.Draw("same")
                mpv = func.GetParameter(1)
                sigma = func.GetParameter(2)
                legendText += f" (MPV={mpv:.0f}, #sigma={sigma:.0f})"
            
            l.AddEntry(h, legendText, "pe")

        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"./plot/{outName}_{suffix}.pdf")

    def plotSeedChargeByCS(self, c, l, datasets: List[DataConfig], suffix: str, plotTitleBase: str):
        valid_datasets = [ds for ds in datasets if ds.histNames.get("seedChargeBySize")]
        if not valid_datasets: return

        for cs in range(1, 7):
            c.Clear(); l.Clear()
            hists = []
            max_y = 0.0
            for ds in valid_datasets:
                baseName = ds.histNames.get("seedChargeBySize")
                targetName = f"{baseName}{cs}"
                h = self.getScaledHist(ds.fileKey, targetName, ds.scale, "", ds.isMerged)
                if h:
                    self.setHistStyle(h, ds.color, ds.markerStyle)

                    rebin_val = ds.rebinValues.get("seedChargeBySize", 1)
                    if rebin_val > 1 and ds.doRebin:
                        h.Rebin(rebin_val)
                    if h.GetMaximum() > 0: h.Scale(1.0 / h.GetMaximum())
                    max_y = max(max_y, h.GetMaximum())
                    hists.append((h, ds.name))

            if not hists: continue
            for i, (h, name) in enumerate(hists):
                opt = "PE" if i == 0 else "same PE"
                if i == 0:
                    h.SetTitle(";charge [e];counts")
                    h.GetXaxis().SetRangeUser(0, 4000)
                    h.GetYaxis().SetRangeUser(0, max_y * 1.2)
                h.Draw(opt)
                l.AddEntry(h, f"{name}, seed", "pe")
            l.Draw()
            self.drawTitle(f"{plotTitleBase} cs{cs}")
            c.SaveAs(f"./plot/seed_charge_cs/clusterSeedCharge_cs{cs}_{suffix}.pdf")

    def drawTitle(self, text):
        t = ROOT.TLatex()
        t.SetTextAlign(12); t.SetTextSize(0.05)
        t.DrawLatexNDC(0.5, 0.85, text)

if __name__ == "__main__":
    plotter = plot_Check()
    try:
        plotter.run_plotCheck()
    except Exception as e:
        log.critical(f"Error: {e}", exc_info=True)
    finally:
        plotter.cleanupResources()