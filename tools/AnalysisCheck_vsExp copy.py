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

# --- ROOTカラーのインポート ---
from ROOT import kBlue, kAzure, kRed, kPink, kGreen, kOrange, kBlack, kCyan

# --- DataConfig ---
@dataclass
class DataConfig:
    """単一のデータソース（実験またはシミュレーション）を定義する"""
    name: str         # 凡例名
    fileKey: str      # ファイルパス or Base
    scale: float      # スケール
    color: int        # 基本色 (ClusterCharge, Comparison用)
    colorAlt: int     # 代替色 (SeedCharge用)
    isMerged: bool    # 実験データか(mergeするか)
    
    # ★ 追加: Rebinを行うかどうか (Sim=False, Exp=Trueを想定)
    doRebin: bool = True 

    histNames: Dict[str, str] = field(default_factory=dict)


# --- ヘルパー関数: get_merged_object ---
def get_merged_object(base_file_path, hist_name):
    # log.info(f"get_merged_object | Merging '{hist_name}'...")
    file_list = glob.glob(f"{base_file_path}_*.root")
    
    if not file_list:
        log.error(f"get_merged_object | No files found for base path: {base_file_path}")
        return None

    merged_hist = None
    
    for file_path in file_list:
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            continue
        
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
            if file:
                file.Close()
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
        hist.SetLineWidth(1) # 線を少し太く

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

        if not hist:
            # log.warning(f"getScaledHist | Not found: {histName}")
            return None

        return self.scalingHistogram(hist, scaleFactor, title)

    def run_plotCheck(self):
        """メイン実行関数"""
        log.info("run_plotCheck | Start")

        # --- データセット定義 ---
        data_entries = {
            # === 実験データ (doRebin=True) ===
            "exp_gap": DataConfig(
                name="Exp GAP",
                fileKey="/home/towa/alice3/hist/sps_check/sps202404_15_gap_10V_SeedThd1000e_NeighborThd200e",
                scale=0.238,
                color=kBlack, colorAlt=13,
                isMerged=True,
                doRebin=True, # ★ 実験データはRebinする
                histNames={
                    "clusterCharge": "AnalysisCE65/CE65_6/cluster/clusterCharge",
                    "seedCharge": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge",
                    "neighborChargeSum": "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum",
                    "clusterSize": "AnalysisCE65/CE65_6/cluster/clusterSize",
                    "seedChargeBySize": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size", 
                }
            ),
            "exp_std": DataConfig(
                name="Exp STD",
                fileKey="/home/towa/alice3/hist/sps_check/sps202404_15_std_10V_SeedThd1000e_NeighborThd200e",
                scale=0.240,
                color=kBlack, colorAlt=13, # 13 is Grey
                isMerged=True,
                doRebin=True, # ★ 実験データはRebinする
                histNames={
                    "clusterCharge": "AnalysisCE65/CE65_6/cluster/clusterCharge",
                    "seedCharge": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge",
                    "neighborChargeSum": "AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum",
                    "clusterSize": "AnalysisCE65/CE65_6/cluster/clusterSize",
                    "seedChargeBySize": "AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size",
                }
            ),
            # === シミュレーションデータ (doRebin=False) ===
            "sim_gap_293k": DataConfig(
                name="Sim GAP (et=1)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_293k_st240e_nt48e_pip_120GeV_masetti_et1.root",
                scale=1000.0,
                color=kRed+1, colorAlt=kPink-2,
                isMerged=False,
                doRebin=False, # ★ SimデータはRebinしない
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_293k": DataConfig(
                name="Sim STD (et=1)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_293k_st240e_nt48e_pip_120GeV_masetti_et1.root",
                scale=1000.0,
                color=kRed+1, colorAlt=kPink-2,
                isMerged=False,
                doRebin=False, # ★ SimデータはRebinしない
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_et0p5": DataConfig(
                name="Sim GAP (et=0.5)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_293k_st240e_nt48e_pip_120GeV_masetti_et0p5.root",
                scale=1000.0,
                color=kAzure+7, colorAlt=kAzure,
                isMerged=False,
                doRebin=False, # ★ SimデータはRebinしない
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_et0p5": DataConfig(
                name="Sim STD (et=0.5)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_293k_st240e_nt48e_pip_120GeV_masetti_et0p5.root",
                scale=1000.0,
                color=kAzure+7, colorAlt=kAzure,
                isMerged=False,
                doRebin=False, # ★ SimデータはRebinしない
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_gap_et2": DataConfig(
                name="Sim GAP (et=2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_293k_st240e_nt48e_pip_120GeV_masetti_.root",
                scale=1000.0,
                color=kGreen+2, colorAlt=kGreen-4,
                isMerged=False,
                doRebin=False, # ★ SimデータはRebinしない
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
            "sim_std_et2": DataConfig(
                name="Sim STD (et=2)",
                fileKey="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_293k_st240e_nt48e_pip_120GeV_masetti_.root",
                scale=1000.0,
                color=kGreen+2, colorAlt=kGreen-4,
                isMerged=False,
                doRebin=False, # ★ SimデータはRebinしない
                histNames={
                    "clusterCharge": "cluster_charge",
                    "seedCharge": "seed_charge",
                    "neighborChargeSum": "cluster_neighbor_charge_sum",
                    "clusterSize": "cluster_size",
                    "seedChargeBySize": "seed_charge_size_",
                }
            ),
        }

        # ★ 変更: 比較リスト定義
        comparisons = [
            # 1. GAPの Exp vs Sim比較
            (["exp_gap", "sim_gap_et0p5", "sim_gap_293k", "sim_gap_et2"], "_gap_default"),
            
            # 2. STDの Exp vs Sim比較
            (["exp_std", "sim_std_et0p5", "sim_std_293k", "sim_std_et2"], "_std_default"),
        ]

        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        legend = ROOT.TLegend(0.5, 0.65, 0.9, 0.8)
        legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.03)

        for keys, suffix in comparisons:
            datasets = [data_entries[k] for k in keys if k in data_entries]
            
            if not datasets:
                log.warning(f"No valid datasets found for suffix {suffix}")
                continue

            log.info(f"Processing comparison set: {keys}")
            
            plot_title = " vs ".join([d.name for d in datasets])
            if len(plot_title) > 50: plot_title = f"Comparison {suffix}"

            self.plotAllCharge(canvas, legend, datasets, suffix, plot_title)
            self.plotClusterSize(canvas, legend, datasets, suffix, plot_title)

            self.plotComparison(canvas, legend, datasets, suffix,
                histKey="seedCharge",
                histLegendSuffix="seed",
                outName="clusterSeedCharge_comp", plotTitle=plot_title,
                rebin=10, xMin=0, xMax=4000)
            
            self.plotComparison(canvas, legend, datasets, suffix,
                histKey="clusterCharge",
                histLegendSuffix="cluster",
                outName="clusterCharge_comp", plotTitle=plot_title,
                rebin=10, xMin=0, xMax=4000)

            self.plotComparison(canvas, legend, datasets, suffix,
                histKey="neighborChargeSum",
                histLegendSuffix="neighbor sum",
                outName="clusterNeighborChargeSum_comp", plotTitle=plot_title,
                rebin=10, xMin=0, xMax=2000)

            self.plotSeedChargeByCS(canvas, legend, datasets, suffix, plot_title)
        
        self.cleanupResources()
        log.info("run_plotCheck | Finished")


    # =================================================================
    # ★ プロット関数群
    # =================================================================

    def plotAllCharge(self, c, l, datasets: List[DataConfig], suffix: str, plotTitle: str):
        """AllCharge (Seed & Cluster) プロット"""
        c.Clear(); l.Clear()

        hists_cl = []
        hists_sd = []
        
        max_y = 0.0

        # 全データセットからヒストグラムを取得してリスト化
        for ds in datasets:
            h_cl = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterCharge"), ds.scale, "", ds.isMerged)
            h_sd = self.getScaledHist(ds.fileKey, ds.histNames.get("seedCharge"), ds.scale, "", ds.isMerged)
            
            if h_cl and h_sd:
                # スタイル設定
                self.setHistStyle(h_cl, ds.color, 20)    # Cluster: 塗りつぶし丸
                self.setHistStyle(h_sd, ds.colorAlt, 24) # Seed: 白抜き丸
                
                # ★ Rebin判定: データセットが許可している場合のみ行う
                if ds.doRebin:
                    h_cl.Rebin(10)
                    h_sd.Rebin(10)
                
                # 正規化 (Normalize to Max 1)
                if h_cl.GetMaximum() > 0: h_cl.Scale(1.0 / h_cl.GetMaximum())
                if h_sd.GetMaximum() > 0: h_sd.Scale(1.0 / h_sd.GetMaximum())

                # 最大Y値の更新
                max_y = max(max_y, h_cl.GetMaximum(), h_sd.GetMaximum())
                
                hists_cl.append((h_cl, ds.name))
                hists_sd.append((h_sd, ds.name))

        if not hists_cl: return

        # 描画ループ
        for i, ((h_cl, name), (h_sd, _)) in enumerate(zip(hists_cl, hists_sd)):
            opt = "PE" if i == 0 else "same PE"
            
            if i == 0:
                h_cl.SetTitle(";charge [e];counts")
                h_cl.GetXaxis().SetRangeUser(0, 4000)
                h_cl.GetYaxis().SetRangeUser(0, max_y * 1.2) # マージン確保
            
            h_cl.Draw(opt)
            h_sd.Draw("same PE")
            
            l.AddEntry(h_sd, f"{name}, seed", "pe")
            l.AddEntry(h_cl, f"{name}, cluster", "pe")

        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"./plot/AllclusterCharge_{suffix}.pdf")

    def plotClusterSize(self, c, l, datasets: List[DataConfig], suffix: str, plotTitle: str):
        """Cluster Size プロット"""
        c.Clear(); l.Clear()

        hists = []
        max_y = 0.0

        for ds in datasets:
            # ClusterSizeは通常Rebinしないので、doRebinフラグはチェックせずそのまま
            h = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterSize"), 1.0, "", ds.isMerged) 
            if h:
                self.setHistStyle(h, ds.color, 20, 0.2) # alpha 0.2
                
                if h.GetEntries() > 0:
                    h.Scale(1.0 / h.GetEntries()) # エントリー数で正規化
                
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
            h.Draw(opt_pe) # 点とヒストグラム両方描画
            
            l.AddEntry(h, f"{name} (mean={h.GetMean():.2f})", "pef")

        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"./plot/clusterSize_{suffix}.pdf")

    def plotComparison(self, c, l, datasets: List[DataConfig], suffix: str,
                       histKey: str, histLegendSuffix: str,
                       outName: str, plotTitle: str,
                       rebin: int, xMin: float, xMax: float):
        """汎用比較プロット"""
        c.Clear(); l.Clear()

        hists = []
        max_y = 0.0

        for ds in datasets:
            histName = ds.histNames.get(histKey)
            if not histName: continue

            h = self.getScaledHist(ds.fileKey, histName, ds.scale, "", ds.isMerged)
            if h:
                self.setHistStyle(h, ds.color, 20)

                # ★ Rebin判定: 引数のrebin > 1 かつ データセットが許可している場合
                if rebin > 1 and ds.doRebin:
                    h.Rebin(rebin)
                
                if h.GetMaximum() > 0:
                    h.Scale(1.0 / h.GetMaximum()) # Normalize to Max 1
                
                max_y = max(max_y, h.GetMaximum())
                hists.append((h, ds.name))

        if not hists: return

        for i, (h, name) in enumerate(hists):
            opt = "PE" if i == 0 else "same PE"
            
            if i == 0:
                h.SetTitle(";charge [e];counts")
                h.GetXaxis().SetRangeUser(xMin, xMax)
                h.GetYaxis().SetRangeUser(0, max_y * 1.2)
            
            h.Draw(opt)
            l.AddEntry(h, f"{name}, {histLegendSuffix}", "pe")

        l.Draw()
        self.drawTitle(plotTitle)
        c.SaveAs(f"./plot/{outName}_{suffix}.pdf")

    def plotSeedChargeByCS(self, c, l, datasets: List[DataConfig], suffix: str, plotTitleBase: str):
        """CSごとのSeed Charge"""
        
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
                    self.setHistStyle(h, ds.color, 24) # Seed style
                    
                    # ★ Rebin判定: データセットが許可している場合のみ
                    if ds.doRebin:
                        h.Rebin(10)

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