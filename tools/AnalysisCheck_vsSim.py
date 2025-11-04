import ROOT
import logging
import os
import glob
from dataclasses import dataclass

# --- ロギング設定 ---
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s | %(levelname)-7s | %(name)-15s | %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger("plot_Check")

# --- ROOTグローバル設定 ---
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
# フィットボックスも非表示に
ROOT.gStyle.SetOptFit(0)

# pixel_pitch = "22p5"

# --- ROOTカラーのインポート ---
from ROOT import kBlue, kAzure, kRed, kPink, kGreen, kOrange, kViolet, kSpring, kCyan, kMagenta, kGray

# --- get_merged_object (変更なし) ---
def get_merged_object(base_file_path, hist_name):
    # (中略: この関数は元のままですが、今回の処理では呼ばれません)
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
        if not isinstance(h_obj, ROOT.TH1):
             log.warning(f"get_merged_object | Object '{hist_name}' in {file_path} is not a TH1")
             f.Close()
             continue
        if merged_hist is None:
            merged_hist = h_obj.Clone(f"{h_obj.GetName()}_merged")
            merged_hist.SetDirectory(0)
        else:
            merged_hist.Add(h_obj)
        f.Close()
    if merged_hist is None:
        log.error(f"get_merged_object | Failed to merge any histograms for '{hist_name}' from base {base_file_path}")
    return merged_hist


class plot_Check:
    def __init__(self):
        """コンストラクタ"""
        self.m_fileCache = {}
        log.info("plot_Check::__init__ | plot_Check object is created.")
        os.makedirs("./plot/seed_charge_cs", exist_ok=True)

        # ★★★ 15個の温度に対応する色のリスト ★★★
        self.colors = [
            kRed + 2, kOrange + 8, kOrange - 2, kSpring + 9, kGreen + 2,
            kCyan + 2, kAzure + 8, kAzure - 4, kBlue + 1, kViolet + 6,
            kViolet - 4, kMagenta + 1, kPink + 6, kPink - 9, kGray + 2
        ]
        # ★★★ 15個の温度に対応するマーカースタイル ★★★
        self.markers = [
            20, 21, 22, 23, 29, 33, 34, 41, 43, 45, 47, 24, 25, 26, 30
        ]

    def __del__(self):
        self.cleanupResources()

    def cleanupResources(self):
        """リソースのクリーンアップ"""
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
        return self.m_fileCache.get(fileName)

    def setHistStyle(self, hist, color, marker, alpha=1.0):
        """ヒストグラムのスタイルを設定"""
        if not hist:
            return
        hist.SetMarkerColor(color)
        hist.SetLineColor(color)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(1.0)
        hist.SetFillColorAlpha(color, alpha)

    def scalingHistogram(self, hist, scale_factor, title):
        """ヒストグラムのX軸を in-place でスケーリングする"""
        if not hist:
            log.warning("scalingHistogram | received a null histogram.")
            return None

        log.info(f"scalingHistogram | Scaling X-axis of histogram: {hist.GetName()} by factor {scale_factor}")

        xAxis = hist.GetXaxis()
        nbins = xAxis.GetNbins()
        xmin = xAxis.GetXmin()
        xmax = xAxis.GetXmax()

        new_xmin = xmin * scale_factor
        new_xmax = xmax * scale_factor

        xAxis.Set(nbins, new_xmin, new_xmax)
        hist.SetTitle(title)

        for i in range(0, nbins + 2):
            content_orig = hist.GetBinContent(i)
            if content_orig < 0:
                hist.SetBinContent(i, 0)

        return hist

    def getScaledHist(self, file_path, histName, scaleFactor, title, isMerged=False):
        """ヒストグラムを取得し、クローンし、X軸をスケーリングする"""
        hist = None
        if isMerged:
            hist = get_merged_object(file_path, histName)
        else:
            f = self.openFile(file_path)
            if f:
                h_obj = f.Get(histName)
                if h_obj:
                    hist = h_obj.Clone(f"{h_obj.GetName()}_clone")
                    hist.SetDirectory(0)
                else:
                    log.warning(f"getScaledHist | Hist '{histName}' not found in file: {file_path}")
            else:
                 log.warning(f"getScaledHist | Could not open file: {file_path}")

        if not hist:
            log.warning(f"getScaledHist | Could not get hist: {histName} from file/base: {file_path}")
            return None

        scaledHist = self.scalingHistogram(hist, scaleFactor, title)
        return scaledHist

    def run_plotCheck(self):
        """★★★ メインの実行関数 (全温度比較) ★★★"""
        log.info("run_plotCheck | Start run_plotCheck (All Temperatures)")

        chip_types = ["std", "gap"]
        #pixel_pitch = "22p5"
        # pitch_configs = [
        #     {"id": "15", "title": "P15"},
        #     {"id": "22p5", "title": "P22.5"}
        # ]

        # ★★★ 比較するすべての温度 ★★★
        temperatures = [273, 293, 313, 333, 353, 373, 393, 413]

        sim_file_base_path = "/home/towa/alice3/hist/sim_modeling/"
        sim_file_template = "analysis_ce65_p22p5_{type}_10v_n0e_{temp}k_st240e_nt48e_pip_120GeV_masetti_.root"

        scale_common = 1000.0

        canvas = ROOT.TCanvas("canvas", "canvas", 900, 700)

        # ★★★ 描画するヒストグラムの設定 (cluster_chargeに "fit": "landau" を追加) ★★★
        hist_configs = [
            {
                "name": "seed_charge",
                "out": "seed_charge_all_temps",
                "title": "Seed Charge",
                "xtitle": "charge [e]", "xmin": 0, "xmax": 4000
            },
            {
                "name": "cluster_charge",
                "out": "cluster_charge_all_temps_WITH_FIT", # ★出力名を変更
                "title": "Cluster Charge (Landau Fit)",       # ★タイトルを変更
                "xtitle": "charge [e]", "xmin": 0, "xmax": 4000,
                "fit": "landau" # ★★★ フィット指定 ★★★
            },
            {
                "name": "cluster_neighbor_charge_sum",
                "out": "neighbor_charge_sum_all_temps",
                "title": "Neighbor Charge Sum",
                "xtitle": "charge [e]", "xmin": 0, "xmax": 2000
            }
        ]

        for chip_type in chip_types:
            # 1. 通常のヒストグラム (全温度比較)
            for config in hist_configs:
                self.plotAllTemperatures(
                    canvas, chip_type, temperatures, config,
                    sim_file_template, sim_file_base_path, scale_common
                )

            # 2. クラスタサイズ別のヒストグラム (全温度比較)
            for cs in range(1, 7):
                self.plotSeedChargeByCS_AllTemps(
                    canvas, chip_type, temperatures, cs,
                    sim_file_template, sim_file_base_path, scale_common
                )

        self.cleanupResources()
        log.info("run_plotCheck | Finished run_plotCheck (All Temperatures)")


    # ★★★ [修正] 全温度を1キャンバスに描画する関数 (Landau fit 修正) ★★★
    def plotAllTemperatures(self, c, chip_type, temperatures, hist_config, file_template, base_path, scale):
        """
        1つのヒストグラムについて、全温度のパターンを1枚のキャンバスに描画する
        """
        log.info(f"Plotting all temperatures for {chip_type.upper()} - {hist_config['name']}")
        c.Clear()

        # 凡例 (Legend)
        legend = ROOT.TLegend(0.60, 0.25, 0.9, 0.9)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.022)

        loaded_items = [] # (hist, temp, fit_function) を格納

        # 1. 全温度のヒストグラムをロード
        for i, temp in enumerate(temperatures):
            file_path = os.path.join(base_path, file_template.format(type=chip_type, temp=temp))
            if not os.path.exists(file_path):
                log.warning(f"File not found, skipping: {file_path}")
                continue

            hist_title = f";{hist_config['xtitle']};Normalized Counts"

            h = self.getScaledHist(file_path, hist_config["name"], scale, hist_title, isMerged=False)
            if not h:
                continue

            color = self.colors[i % len(self.colors)]
            marker = self.markers[i % len(self.markers)]
            self.setHistStyle(h, color, marker)

            fit_clone = None
            norm_factor = h.GetMaximum()

            # ★★★ Landau Fit 処理 (修正) ★★★
            if hist_config.get("fit") == "landau":
                if norm_factor > 0:
                    # フィットを実行 (正規化前)
                    # ★★★ "NQ" (保存しない) を "Q0" (描画しないが保存する) に変更 ★★★
                    h.Fit("landau", "Q0")
                    fit_func = h.GetFunction("landau") # これで関数を取得できる

                    if fit_func:
                        log.info(f"Fit successful for {chip_type} {temp}K.")
                        # スタイル設定用にクローン
                        fit_clone = fit_func.Clone(f"{h.GetName()}_fit_clone")
                        fit_clone.SetLineColor(color)
                        fit_clone.SetLineStyle(2) # 破線
                        fit_clone.SetLineWidth(2)

                        # ★フィット関数も正規化
                        # (Landauの Constant は Par[0])
                        fit_clone.SetParameter(0, fit_clone.GetParameter(0) / norm_factor)
                    else:
                        log.warning(f"Fit failed for {chip_type} {temp}K, GetFunction returned None.")
                else:
                    log.warning(f"Histogram {hist_config['name']} for temp {temp} is empty, skipping fit.")


            # ★ 形状比較のため、各ヒストグラムを最大値で正規化 (フィット後)
            if norm_factor > 0:
                h.Scale(1.0 / norm_factor)
            else:
                log.warning(f"Histogram {hist_config['name']} for temp {temp} is empty.")
                continue

            h.GetYaxis().SetRangeUser(0, 1.2) # Y軸を固定
            h.GetXaxis().SetRangeUser(hist_config["xmin"], hist_config["xmax"])

            loaded_items.append((h, temp, fit_clone)) # ★フィット関数も追加

        if not loaded_items:
            log.error(f"No histograms loaded for {chip_type} - {hist_config['name']}")
            return

        # 2. 描画
        # 最初にダミーを描画して軸を設定
        h_first, _, _ = loaded_items[0]
        h_first.Draw("PE")

        for i, (h, temp, fit_func) in enumerate(loaded_items):
            draw_opt = "SAME PE" # 常に重ね描き
            if i == 0:
                draw_opt = "PE" # 1枚目だけ

            h.Draw(draw_opt)
            legend.AddEntry(h, f"{chip_type.upper()} {temp}K", "pe")

            # ★フィット関数を描画
            if fit_func:
                fit_func.Draw("SAME L")
                mpv = fit_func.GetParameter(1) # Par[1] が MPV
                # ★凡例にフィット情報を追加
                legend.AddEntry(fit_func, f"  MPV: {mpv:.1f} e", "l")


        legend.Draw()

        # 3. タイトル描画
        plot_title = f"SQ P22.5 {chip_type.upper()} 10V - {hist_config['title']}"
        latexTitle = ROOT.TLatex()
        latexTitle.SetNDC()
        latexTitle.SetTextAlign(22)
        latexTitle.SetTextSize(0.04)
        latexTitle.DrawLatex(0.5, 0.95, plot_title)

        # 4. 保存
        output_filename = f"./plot/{hist_config['out']}_{chip_type}.pdf"
        c.SaveAs(output_filename)
        log.info(f"Saved plot: {output_filename}")

    # ★★★ [変更なし] CS別に全温度を1キャンバスに描画する関数 ★★★
    def plotSeedChargeByCS_AllTemps(self, c, chip_type, temperatures, cs, file_template, base_path, scale):
        """
        CS（クラスタサイズ）別に、全温度のSeed Chargeを1枚のキャンバスに描画する
        """
        hist_name = f"seed_charge_size_{cs}"
        log.info(f"Plotting all temperatures for {chip_type.upper()} - {hist_name}")
        c.Clear()

        # 凡例のテキストサイズを、色/マーカーの数に合わせて調整
        legend_text_size = 0.025
        if len(temperatures) > 12:
             legend_text_size = 0.022

        legend = ROOT.TLegend(0.65, 0.45, 0.9, 0.9)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(legend_text_size)

        loaded_hists = []

        # 1. ロード
        for i, temp in enumerate(temperatures):
            file_path = os.path.join(base_path, file_template.format(type=chip_type, temp=temp))
            if not os.path.exists(file_path):
                log.warning(f"File not found, skipping: {file_path}")
                continue

            hist_title = ";charge [e];Normalized Counts"

            h = self.getScaledHist(file_path, hist_name, scale, hist_title, isMerged=False)
            if not h:
                continue

            if h.GetMaximum() > 0:
                h.Scale(1.0 / h.GetMaximum())
            else:
                log.warning(f"Histogram {hist_name} for temp {temp} is empty.")
                continue

            color = self.colors[i % len(self.colors)]
            marker = self.markers[i % len(self.markers)]
            self.setHistStyle(h, color, marker)

            h.GetYaxis().SetRangeUser(0, 1.2)
            h.GetXaxis().SetRangeUser(0, 4000) # CS別は 0-4000 で固定

            loaded_hists.append((h, temp))

        if not loaded_hists:
            log.error(f"No histograms loaded for {chip_type} - {hist_name}")
            return

        # 2. 描画
        h_first, _ = loaded_hists[0]
        h_first.Draw("PE") # 軸の設定用

        for i, (h, temp) in enumerate(loaded_hists):
            draw_opt = "SAME PE"
            if i == 0:
                draw_opt = "PE"
            h.Draw(draw_opt)
            legend.AddEntry(h, f"{chip_type.upper()} {temp}K", "pe")

        legend.Draw()

        # 3. タイトル
        plot_title = f"SQ P22.5 {chip_type.upper()} 10V - Seed Charge (CS={cs})"
        latexTitle = ROOT.TLatex()
        latexTitle.SetNDC()
        latexTitle.SetTextAlign(22)
        latexTitle.SetTextSize(0.04)
        latexTitle.DrawLatex(0.5, 0.95, plot_title)

        # 4. 保存
        output_filename = f"./plot/seed_charge_cs/seed_charge_cs{cs}_all_temps_{chip_type}.pdf"
        c.SaveAs(output_filename)
        log.info(f"Saved plot: {output_filename}")


# --- スクリプトの実行 ---
if __name__ == "__main__":
    plotter = plot_Check()
    try:
        plotter.run_plotCheck()

    except Exception as e:
        log.critical(f"An unhandled exception occurred: {e}", exc_info=True)
    finally:
        plotter.cleanupResources()
        log.info("Main script finished.")