import ROOT
import logging
import os
import glob
from dataclasses import dataclass

# --- ロギング設定 (C++のLOG_STATUS/LOG_ERRORの代わり) ---
# ログレベルをINFOに設定
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s | %(levelname)-7s | %(name)-15s | %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger("plot_Check")

# --- ROOTグローバル設定 ---
ROOT.gROOT.SetBatch(True) # バッチモードで実行 (GUIウィンドウを表示しない)

# ★★★ 修正点 ★★★
# 統計ボックスをデフォルトで非表示にする
ROOT.gStyle.SetOptStat(0) 

# --- ROOTカラーのインポート (kBlue, kRedなど) ---
from ROOT import kBlue, kAzure, kRed, kPink

# --- C++のDatasetConfig構造体の代わり ---
@dataclass
class DatasetConfig:
    name: str
    suffix: str
    expFileBase: str
    expScale: float
    simFile: str
    simScale: float

# --- plot_BeamTest::get_merged_object のPython版 (推測に基づく実装) ---
def get_merged_object(base_file_path, hist_name):
    """
    base_file_path_*.root に一致する全ファイルから
    hist_name のヒストグラムをマージして返す。
    """
    log.info(f"get_merged_object | Merging '{hist_name}' from files matching '{base_file_path}_*.root'")
    file_list = glob.glob(f"{base_file_path}_*.root")
    
    if not file_list:
        log.error(f"get_merged_object | No files found for base path: {base_file_path}")
        return None

    merged_hist = None
    
    for file_path in file_list:
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            log.warning(f"get_merged_object | Failed to open file: {file_path}")
            continue
        
        h_obj = f.Get(hist_name)
        if not h_obj:
            log.warning(f"get_merged_object | Hist '{hist_name}' not found in {file_path}")
            f.Close()
            continue
        
        # オブジェクトがTH1から派生しているか確認
        if not isinstance(h_obj, ROOT.TH1):
             log.warning(f"get_merged_object | Object '{hist_name}' in {file_path} is not a TH1")
             f.Close()
             continue

        if merged_hist is None:
            # 最初のヒストグラムをクローン
            merged_hist = h_obj.Clone(f"{h_obj.GetName()}_merged")
            merged_hist.SetDirectory(0) # ファイルから切り離す
        else:
            # 2つ目以降を足し合わせる
            merged_hist.Add(h_obj)
        
        f.Close() # ファイルを閉じる
    
    if merged_hist is None:
        log.error(f"get_merged_object | Failed to merge any histograms for '{hist_name}' from base {base_file_path}")

    return merged_hist


class plot_Check:
    def __init__(self):
        """コンストラクタ"""
        self.m_fileCache = {}  # C++のm_fileCache (std::map) の代わり
        log.info("plot_Check::__init__ | plot_Check object is created.")
        
        # プロット保存用ディレクトリを作成
        os.makedirs("./plot/seed_charge_cs", exist_ok=True)


    def __del__(self):
        """デストラクタ (Pythonでは__del__は非推奨だが、念のため)"""
        self.cleanupResources()

    def cleanupResources(self):
        """リソースのクリーンアップ (C++のデストラクタ内の処理)"""
        for file_name, file in self.m_fileCache.items():
            if file:
                log.info(f"cleanupResources | Closing file: {file_name}")
                file.Close()
        self.m_fileCache.clear()
        log.info("cleanupResources | Cleaned up file cache.")

    def openFile(self, fileName):
        """ファイルを開き、キャッシュに保存する"""
        if fileName not in self.m_fileCache:
            file = ROOT.TFile.Open(fileName)
            if not file or file.IsZombie():
                log.error(f"openFile | Failed to open file: {fileName}")
                self.m_fileCache[fileName] = None
            else:
                log.info(f"openFile | Opened and cached file: {fileName}")
                self.m_fileCache[fileName] = file
        return self.m_fileCache.get(fileName) # .get() で安全に取得

    def setHistStyle(self, hist, color, marker, alpha=1.0):
        """ヒストグラムのスタイルを設定"""
        if not hist: 
            return
        hist.SetMarkerColor(color)
        hist.SetLineColor(color)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(0.8)
        hist.SetFillColorAlpha(color, alpha)

    def scalingHistogram(self, hist, scale_factor, title):
        """
        ヒストグラムのX軸を **in-place (元のオブジェクトを直接)** でスケーリングする
        C++版のロジックと同一。
        """
        if not hist:
            log.warning("scalingHistogram | received a null histogram.")
            return None
        
        log.info(f"scalingHistogram | Scaling histogram: {hist.GetName()} by factor {scale_factor}")

        xAxis = hist.GetXaxis()
        nbins = xAxis.GetNbins()
        xmin = xAxis.GetXmin()
        xmax = xAxis.GetXmax()

        # 新しいX軸の範囲を計算
        new_xmin = xmin * scale_factor
        new_xmax = xmax * scale_factor

        # ヒストグラムの軸を直接設定 (in-place)
        xAxis.Set(nbins, new_xmin, new_xmax)
        
        # タイトルと軸ラベルを更新
        hist.SetTitle(title)

        # ビンの中身をチェック (0からnbins+1まで)
        for i in range(0, nbins + 2):
            content_orig = hist.GetBinContent(i)
            if content_orig < 0:
                hist.SetBinContent(i, 0)
                
        return hist # 元のヒストグラムへの参照を返す

    def getScaledHist(self, file_path, histName, scaleFactor, title, isMerged=False):
        """ヒストグラムを取得し、クローンし、スケーリングする"""
        hist = None
        if isMerged:
            # 外部ヘルパー関数 (plot_BeamTest::get_merged_objectの代わり) を呼ぶ
            hist = get_merged_object(file_path, histName)
        else:
            f = self.openFile(file_path)
            if f:
                h_obj = f.Get(histName)
                if h_obj:
                    # クローンを作成 (重要)
                    hist = h_obj.Clone(f"{h_obj.GetName()}_clone")
                    hist.SetDirectory(0) # 元のファイルから切り離す
                else:
                    log.warning(f"getScaledHist | Hist '{histName}' not found in file: {file_path}")
            else:
                 log.warning(f"getScaledHist | Could not open file: {file_path}")

        if not hist:
            log.warning(f"getScaledHist | Could not get hist: {histName} from file/base: {file_path}")
            return None

        # scalingHistogram は渡されたヒストグラム(hist)を直接変更する
        scaledHist = self.scalingHistogram(hist, scaleFactor, title)
        
        # scaledHist と hist は同じオブジェクトを指している
        return scaledHist

    def run_plotCheck(self):
        """メインの実行関数"""
        log.info("run_plotCheck | Start run_plotCheck")

        datasets = [
            DatasetConfig(
                name="SQ P22.5 GAP 10V", 
                suffix="_gap",
                expFileBase="/home/towa/alice3/hist/sps_check/sps202404_22p5_gap_10V_SeedThd1000e_NeighborThd200e", 
                expScale=0.238,
                simFile="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p22p5_gap_10v_n0e_323k_st240e_nt48e_pip_120GeV_masetti_.root", 
                simScale=1000.0
            ),
            DatasetConfig(
                name="SQ P22.5 STD 10V", 
                suffix="_std",
                expFileBase="/home/towa/alice3/hist/sps_check/sps202404_22p5_std_10V_SeedThd1000e_NeighborThd200e", 
                expScale=0.240,
                simFile="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p22p5_std_10v_n0e_323k_st240e_nt48e_pip_120GeV_masetti_.root", 
                simScale=1000.0
            )
        ]
        # datasets = [
        #     DatasetConfig(
        #         name="SQ P15 GAP 10V", 
        #         suffix="_gap",
        #         expFileBase="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_293_st240e_nt48e_pip_120GeV_masetti", 
        #         expScale=1000.0,
        #         simFile="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_gap_10v_n0e_373_st240e_nt48e_pip_120GeV_masetti_.root", 
        #         simScale=1000.0
        #     ),
        #     DatasetConfig(
        #         name="SQ P15 STD 10V", 
        #         suffix="_std",
        #         expFileBase="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_293_st240e_nt48e_pip_120GeV_masetti", 
        #         expScale=1000.0,
        #         simFile="/home/towa/alice3/hist/sim_modeling/analysis_ce65_p15_std_10v_n0e_373_st240e_nt48e_pip_120GeV_masetti_.root", 
        #         simScale=1000.0
        #     )
        # ]

        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        legend = ROOT.TLegend(0.5, 0.65, 0.9, 0.8)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.03)

        for ds in datasets:
            log.info(f"Processing dataset: {ds.name}")

            self.plotAllCharge(canvas, legend, ds)
            self.plotClusterSize(canvas, legend, ds)

            self.plotComparison(canvas, legend, ds,
                expHistName="AnalysisCE65/CE65_6/cluster/clusterSeedCharge", simHistName="seed_charge",
                expLegend="exp, seed", simLegend="sim, seed",
                outName="clusterSeedCharge_sim_exp", plotTitle=ds.name,
                rebin=10, xMin=0, xMax=4000)
            
            self.plotComparison(canvas, legend, ds,
                expHistName="AnalysisCE65/CE65_6/cluster/clusterCharge", simHistName="cluster_charge",
                expLegend="exp, cluster", simLegend="sim, cluster",
                outName="clusterCharge_sim_exp", plotTitle=ds.name,
                rebin=10, xMin=0, xMax=4000)

            self.plotComparison(canvas, legend, ds,
                expHistName="AnalysisCE65/CE65_6/cluster/clusterNeighborsChargeSum", simHistName="cluster_neighbor_charge_sum",
                expLegend="exp, neighbor sum", simLegend="sim, neighbor sum",
                outName="clusterNeighborChargeSum_sim_exp", plotTitle=ds.name,
                rebin=10, xMin=0, xMax=2000)

            # self.plotComparison(canvas, legend, ds,
            #     expHistName="seed_charge", simHistName="seed_charge",
            #     expLegend="293K, seed", simLegend="373K, seed",
            #     outName="clusterSeedCharge_sim_exp", plotTitle=ds.name,
            #     rebin=1, xMin=0, xMax=4000)
            
            # self.plotComparison(canvas, legend, ds,
            #     expHistName="cluster_charge", simHistName="cluster_charge",
            #     expLegend="293K, cluster", simLegend="373K, cluster",
            #     outName="clusterCharge_sim_exp", plotTitle=ds.name,
            #     rebin=1, xMin=0, xMax=4000)

            # self.plotComparison(canvas, legend, ds,
            #     expHistName="cluster_neighbor_charge_sum", simHistName="cluster_neighbor_charge_sum",
            #     expLegend="293K, neighbor sum", simLegend="373K, neighbor sum",
            #     outName="clusterNeighborChargeSum_sim_exp", plotTitle=ds.name,
            #     rebin=1, xMin=0, xMax=2000)

            self.plotSeedChargeByCS(canvas, legend, ds)
        
        # C++の `delete canvas; delete legend;` はPythonでは不要。
        # スコープを抜ければガベージコレクタが処理する。
        # 最後にファイルキャッシュを明示的にクリーンアップする
        self.cleanupResources()

        log.info("run_plotCheck | Finished run_plotCheck")

    def plotAllCharge(self, c, l, ds):
        """AllChargeプロットを作成"""
        log.info(f"plotAllCharge | Plotting all charge for {ds.name}")
        c.Clear()
        l.Clear()

        h_exp_cl = self.getScaledHist(ds.expFileBase, "AnalysisCE65/CE65_6/cluster/clusterCharge", ds.expScale, ";charge [e];counts", isMerged=True)
        h_exp_sd = self.getScaledHist(ds.expFileBase, "AnalysisCE65/CE65_6/cluster/clusterSeedCharge", ds.expScale, ";charge [e];counts", isMerged=True)
        h_sim_cl = self.getScaledHist(ds.simFile, "cluster_charge", ds.simScale, ";charge [e];counts", isMerged=False)
        h_sim_sd = self.getScaledHist(ds.simFile, "seed_charge", ds.simScale, ";charge [e];counts", isMerged=False)

        # h_exp_cl = self.getScaledHist(ds.expFileBase, "clsuter_charge", ds.expScale, ";charge [e];counts", isMerged=True)
        # h_exp_sd = self.getScaledHist(ds.expFileBase, "seed_charge", ds.expScale, ";charge [e];counts", isMerged=True)
        # h_sim_cl = self.getScaledHist(ds.simFile, "cluster_charge", ds.simScale, ";charge [e];counts", isMerged=False)
        # h_sim_sd = self.getScaledHist(ds.simFile, "seed_charge", ds.simScale, ";charge [e];counts", isMerged=False)

        # 4つのヒストグラムがすべて存在するかチェック
        if not all([h_exp_cl, h_exp_sd, h_sim_cl, h_sim_sd]):
            log.error("plotAllCharge | Missing one or more histograms for plotAllCharge. Skipping.")
            return

        self.setHistStyle(h_exp_cl, kBlue, 20)
        self.setHistStyle(h_exp_sd, kAzure - 3, 24)
        self.setHistStyle(h_sim_cl, kRed + 1, 20)
        self.setHistStyle(h_sim_sd, kPink - 2, 24)

        h_exp_cl.Rebin(10)
        h_exp_sd.Rebin(10)

        # 正規化
        if h_exp_cl.GetMaximum() > 0: h_exp_cl.Scale(1.0 / h_exp_cl.GetMaximum())
        if h_exp_sd.GetMaximum() > 0: h_exp_sd.Scale(1.0 / h_exp_sd.GetMaximum())
        if h_sim_cl.GetMaximum() > 0: h_sim_cl.Scale(1.0 / h_sim_cl.GetMaximum())
        if h_sim_sd.GetMaximum() > 0: h_sim_sd.Scale(1.0 / h_sim_sd.GetMaximum())

        h_sim_cl.GetXaxis().SetRangeUser(0, 4000)
        h_sim_cl.GetYaxis().SetRangeUser(0, 1.1)

        h_sim_cl.Draw("PE")
        h_exp_cl.Draw("samePE")
        h_sim_sd.Draw("samePE")
        h_exp_sd.Draw("samePE")

        l.AddEntry(h_exp_sd, "exp, seed", "pe")
        l.AddEntry(h_exp_cl, "exp, cluster", "pe")
        l.AddEntry(h_sim_sd, "sim, seed", "pe")
        l.AddEntry(h_sim_cl, "sim, cluster", "pe")
        l.Draw()

        title = ROOT.TLatex()
        title.SetTextAlign(12)
        title.SetTextSize(0.05)
        title.DrawLatexNDC(0.6, 0.85, ds.name)

        c.SaveAs(f"./plot/AllclusterCharge_sim_exp{ds.suffix}.pdf")
        # Pythonではローカル変数のヒストグラムの 'delete' は不要

    def plotClusterSize(self, c, l, ds):
        """ClusterSizeプロットを作成"""
        log.info(f"plotClusterSize | Plotting ClusterSize for {ds.name}")
        c.Clear()
        l.Clear()

        h_exp_size_orig = get_merged_object(ds.expFileBase, "AnalysisCE65/CE65_6/cluster/clusterSize")
        # h_exp_size_orig = get_merged_object(ds.expFileBase, "cluster_size")
        

        h_sim_size_orig = None
        f_sim = self.openFile(ds.simFile)
        if f_sim:
            h_sim_size_orig = f_sim.Get("cluster_size")

        if not h_exp_size_orig or not h_sim_size_orig:
            log.error("plotClusterSize | Missing histograms for plotClusterSize. Skipping.")
            return
        
        # C++コードと同様に、取得後にクローンする
        h_exp_size = h_exp_size_orig.Clone(f"{h_exp_size_orig.GetName()}_clone")
        h_exp_size.SetDirectory(0)
        h_sim_size = h_sim_size_orig.Clone(f"{h_sim_size_orig.GetName()}_clone")
        h_sim_size.SetDirectory(0)

        self.setHistStyle(h_exp_size, kBlue, 20, 0.2)
        self.setHistStyle(h_sim_size, kRed + 1, 20, 0.2)

        # エントリー数で正規化
        if h_exp_size.GetEntries() > 0: h_exp_size.Scale(1.0 / h_exp_size.GetEntries())
        if h_sim_size.GetEntries() > 0: h_sim_size.Scale(1.0 / h_sim_size.GetEntries())
        
        h_sim_size.GetXaxis().SetRangeUser(1, 10)
        h_sim_size.GetYaxis().SetRangeUser(0, 1)

        h_sim_size.Draw("HIST")
        h_exp_size.Draw("same HIST")
        h_sim_size.Draw("same PE")
        h_exp_size.Draw("samePE")

        l.AddEntry(h_exp_size, "exp", "pef")
        l.AddEntry(h_exp_size, f"mean = {h_exp_size.GetMean():.2f}", "")
        l.AddEntry(h_sim_size, "sim", "pef")
        l.AddEntry(h_sim_size, f"mean = {h_sim_size.GetMean():.2f}", "")
        l.Draw()

        title = ROOT.TLatex()
        title.SetTextAlign(12)
        title.SetTextSize(0.05)
        title.DrawLatexNDC(0.6, 0.85, ds.name)
        
        c.SaveAs(f"./plot/clusterSize_sim_exp{ds.suffix}.pdf")

    def plotComparison(self, c, l, ds,
                       expHistName, simHistName,
                       expLegend, simLegend,
                       outName, plotTitle,
                       rebin, xMin, xMax, normalizeToMax=True):
        """汎用的な比較プロット関数"""
        log.info(f"plotComparison | Plotting Comparison {outName} for {ds.name}")
        c.Clear()
        l.Clear()

        title = ";charge [e];counts"
        h_exp = self.getScaledHist(ds.expFileBase, expHistName, ds.expScale, title, isMerged=True)
        h_sim = self.getScaledHist(ds.simFile, simHistName, ds.simScale, title, isMerged=False)

        if not h_exp or not h_sim:
            log.error(f"plotComparison | Missing histograms for {outName}. Skipping.")
            return

        self.setHistStyle(h_exp, kBlue, 20)
        self.setHistStyle(h_sim, kRed + 1, 20)

        if rebin > 1:
            h_exp.Rebin(rebin)
            # C++コードはh_expのみRebinしているため、それに倣う
        
        if normalizeToMax:
            if h_exp.GetMaximum() > 0: h_exp.Scale(1.0 / h_exp.GetMaximum())
            if h_sim.GetMaximum() > 0: h_sim.Scale(1.0 / h_sim.GetMaximum())

        h_sim.GetXaxis().SetRangeUser(xMin, xMax)
        
        # Y軸の範囲を自動調整 (もしくは固定)
        if normalizeToMax:
            h_sim.GetYaxis().SetRangeUser(0, 1.1)
        else:
            maxY = max(h_sim.GetMaximum(), h_exp.GetMaximum())
            h_sim.GetYaxis().SetRangeUser(0, maxY * 1.1)

        h_sim.Draw("PE")
        h_exp.Draw("samePE")

        l.AddEntry(h_exp, expLegend, "pe")
        l.AddEntry(h_sim, simLegend, "pe")
        l.Draw()
        
        latexTitle = ROOT.TLatex()
        latexTitle.SetTextAlign(12)
        latexTitle.SetTextSize(0.05)
        latexTitle.DrawLatexNDC(0.6, 0.85, plotTitle)

        c.SaveAs(f"./plot/{outName}{ds.suffix}.pdf")

    def plotSeedChargeByCS(self, c, l, ds):
        """クラスタサイズ別のSeed Chargeプロット"""
        log.info(f"plotSeedChargeByCS | Plotting SeedChargeByCS for {ds.name}")

        for cs in range(1, 7): # 1から6まで
            c.Clear()
            l.Clear()

            expHistName = f"AnalysisCE65/CE65_6/cluster/clusterSeedCharge_size{cs}"
            # C++コードのForm("..._%d", cs) は "..._1" となる
            simHistName = f"seed_charge_size_{cs}" 
            plotTitle = f"{ds.name} cs{cs}"
            outName = f"seed_charge_cs/clusterSeedCharge_cs{cs}_sim_exp"

            title = ";charge [e];counts"
            h_exp = self.getScaledHist(ds.expFileBase, expHistName, ds.expScale, title, isMerged=True)
            h_sim = self.getScaledHist(ds.simFile, simHistName, ds.simScale, title, isMerged=False)

            if not h_exp or not h_sim:
                log.warning(f"plotSeedChargeByCS | Missing histograms for {outName}. Skipping cs={cs}.")
                continue # 次のクラスタサイズへ

            self.setHistStyle(h_exp, kBlue, 24)
            self.setHistStyle(h_sim, kRed, 24)
            
            h_exp.Rebin(10)

            # 正規化
            if h_sim.GetMaximum() > 0: h_sim.Scale(1.0 / h_sim.GetMaximum())
            if h_exp.GetMaximum() > 0: h_exp.Scale(1.0 / h_exp.GetMaximum())

            h_sim.GetXaxis().SetRangeUser(0, 4000)
            h_sim.GetYaxis().SetRangeUser(0, 1.1) # Y軸を固定
            
            h_sim.Draw("PE")
            h_exp.Draw("samePE")

            l.AddEntry(h_exp, "exp, seed", "pe")
            l.AddEntry(h_sim, "sim, seed", "pe")
            l.Draw()

            latexTitle = ROOT.TLatex()
            latexTitle.SetTextAlign(12)
            latexTitle.SetTextSize(0.05)
            latexTitle.DrawLatexNDC(0.5, 0.85, plotTitle)

            c.SaveAs(f"./plot/{outName}{ds.suffix}.pdf")

    def plotDepositedCharge(self, c, l):
        """DepositedChargeプロット (run_plotCheckからは呼ばれない)"""
        log.info("plotDepositedCharge | Plotting DepositedCharge")
        c.Clear()
        l.Clear()

        file_path = "/home/towa/alice3/hist/test/CE65_sq_p22p5_gap_10v_pip_120GeV_default.root"
        
        h_dep = self.getScaledHist(file_path, "DepositionGeant4/deposited_charge_CE65", 1000, ";charge [e];counts", isMerged=False)
        h_cls = self.getScaledHist(file_path, "DetectorHistogrammer/CE65/charge/cluster_charge", 1000, ";charge [e];counts", isMerged=False)
        
        if not h_dep or not h_cls:
            log.error("plotDepositedCharge | Missing histograms for plotDepositedCharge. Skipping.")
            return

        self.setHistStyle(h_dep, kBlue, 20)
        self.setHistStyle(h_cls, kRed + 1, 20)
        
        h_dep.GetXaxis().SetRangeUser(0, 20000)
        
        # Y軸の範囲を自動設定
        maxY = max(h_dep.GetMaximum(), h_cls.GetMaximum())
        h_dep.GetYaxis().SetRangeUser(0, maxY * 1.1)
        
        h_dep.Draw("PE")
        h_cls.Draw("samePE")

        l.AddEntry(h_dep, "deposited charge", "pe")
        l.AddEntry(h_dep, f"mean = {h_dep.GetMean():.2f}", "")
        l.AddEntry(h_cls, "cluster charge", "pe")
        l.AddEntry(h_cls, f"mean = {h_cls.GetMean():.2f}", "")
        l.Draw()

        c.SaveAs("./plot/depositedCharge.pdf")

# --- スクリプトの実行 ---
if __name__ == "__main__":
    plotter = plot_Check()
    try:
        plotter.run_plotCheck()
        
    except Exception as e:
        log.critical(f"An unhandled exception occurred: {e}", exc_info=True)
    finally:
        # スクリプト終了時に明示的にリソースを解放する
        plotter.cleanupResources()
        log.info("Main script finished.")