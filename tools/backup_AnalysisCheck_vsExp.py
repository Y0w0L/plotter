import ROOT
import logging
import os
import glob
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import json
import math

# --- Langau関数の定義 (C++マクロの埋め込み) ---
# パラメータ:
# [0] Width (Scale) of Landau
# [1] MPV (Most Probable Value)
# [2] Total Area (Integral)
# [3] Sigma of Gaussian (Noise)
langau_code = """
#include <TMath.h>
#include <TF1.h>

Double_t langaufun(Double_t *x, Double_t *par) {
   //Fit parameters:
   //par[0]=Width (scale) of Landau (MPV)
   //par[1]=Most Probable (MPV)
   //par[2]=Total Area (Integral - not the height)
   //par[3]=Width (Sigma) of convoluted Gaussian function

   // Numeric constants
   Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
   Double_t mpshift  = -0.22278298;       // Landau maximum location

   // Control constants
   Double_t np = 100.0;      // number of convolution steps
   Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

   // Variables
   Double_t xx;
   Double_t mpc;
   Double_t fland;
   Double_t sum = 0.0;
   Double_t xlow,xupp;
   Double_t step;
   Double_t i;

   // MP shift correction
   mpc = par[1] - mpshift * par[0];

   // Range of convolution integral
   xlow = x[0] - sc * par[3];
   xupp = x[0] + sc * par[3];

   step = (xupp-xlow) / np;

   // Convolution Loop
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
    rebinType: str = "default"

def parse_root_color(color_str: str) -> int:
    """
    "kRed+1" や "kAzure-4" のような文字列を解析して ROOTの色コード(int)に変換する。
    evalを使わず文字列操作で行うため安全かつ確実。
    """
    if isinstance(color_str, int):
        return color_str
    
    # 空白を削除
    color_str = color_str.replace(" ", "")
    
    try:
        # 加算の場合 (例: "kRed+1")
        if "+" in color_str:
            parts = color_str.split("+")
            base_color_name = parts[0]
            offset = int(parts[1])
            return getattr(ROOT, base_color_name) + offset
            
        # 減算の場合 (例: "kRed-4")
        elif "-" in color_str:
            parts = color_str.split("-")
            base_color_name = parts[0]
            offset = int(parts[1])
            return getattr(ROOT, base_color_name) - offset
            
        # そのままの場合 (例: "kBlack")
        else:
            return getattr(ROOT, color_str)
            
    except AttributeError:
        log.error(f"Color name '{color_str}' not found in ROOT module.")
        return ROOT.kBlack # デフォルト
    except ValueError:
        log.error(f"Invalid color format '{color_str}'. Expected format like 'kRed+1'.")
        return ROOT.kBlack
    except Exception as e:
        log.error(f"Failed to parse color string '{color_str}': {e}")
        return ROOT.kBlack

def load_data_config(json_path: str) -> Dict[str, DataConfig]:
    """
    JSONファイルからDataConfigの辞書を生成する
    """
    if not os.path.exists(json_path):
        log.error(f"Config file not found: {json_path}")
        return {}
    
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    configs = {}
    for key, val in data.items():
        if "color" in val:
            val["color"] = parse_root_color(val["color"])
        if "colorAlt" in val:
            val["colorAlt"] = parse_root_color(val["colorAlt"])
        
        try:
            configs[key] = DataConfig(**val)
        except TypeError as e:
            log.error(f"Error creating DataConfig for '{key}': {e}")
        
    return configs


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
    def __init__(self, config_file="simulation_data.json"):
        self.m_fileCache = {} 
        self.config_file = config_file
        log.info("plot_Check::__init__ | plot_Check object is created.")
        
        # --- Change: Create sub-directories for organized output ---
        self.output_dirs = {
            "1D": "./plot/1D_Comparison",
            "2D": "./plot/2D_Maps",
            "Path": "./plot/PathAnalysis",
            "SeedCS": "./plot/1D_Comparison/seed_charge_cs"
        }
        for d in self.output_dirs.values():
            os.makedirs(d, exist_ok=True)
    # def __init__(self, config_file="simulation_data.json"):
    #     self.m_fileCache = {} 
    #     self.config_file = config_file
    #     log.info("plot_Check::__init__ | plot_Check object is created.")
    #     os.makedirs("./plot/seed_charge_cs", exist_ok=True)

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

    def setHistStyle(self, hist, color, marker, alpha=0.15):
        if not hist: return
        hist.SetMarkerColor(color)
        hist.SetLineColor(color)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(0.8)
        hist.SetFillColorAlpha(color, alpha)
        hist.SetLineWidth(2)

    def setAxisStyle(self, hist):
        """軸のラベルサイズとタイトルサイズを一括設定"""
        if not hist: return
        
        # サイズ設定 (適宜数値を調整してください)
        label_size = 0.04   # 数値の文字サイズ
        title_size = 0.05   # タイトルの文字サイズ
        
        # X軸設定
        xaxis = hist.GetXaxis()
        xaxis.SetLabelSize(label_size)
        xaxis.SetTitleSize(title_size)
        xaxis.SetTitleOffset(0.9) # タイトルと軸の距離
        
        # Y軸設定
        yaxis = hist.GetYaxis()
        yaxis.SetLabelSize(label_size)
        yaxis.SetTitleSize(title_size)
        yaxis.SetTitleOffset(1) # Y軸は数字と被らないよう少し離す

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

        is_neighbor_hist = ("neighbor" in histName.lower())
        is_exp_data = ("sps" in fileKey)

        if is_neighbor_hist and is_exp_data:
            # カットする閾値 (ADC値)
            threshold_adc = 100.0
            
            nbins = hist.GetNbinsX()
            # アンダーフロー(0)からオーバーフロー(nbins+1)まで走査
            for i in range(0, nbins + 2):
                bin_center = hist.GetBinCenter(i)
                
                # ビンの中心値が閾値未満なら内容を消去
                if bin_center < threshold_adc:
                    hist.SetBinContent(i, 0)
                    hist.SetBinError(i, 0)
            
            log.info(f"getScaledHist | Applied cut < {threshold_adc} (ADC) for {histName} before scaling.")
        # --------------------------------------------------------------------

        # --- ★★★ Neighbor Charge Sum の 0 エントリー除外処理 ★★★ ---
        # cluster_neighbor_charge_sum ヒストグラム、またはそれに類する名前の場合
        #if "clusterNeighborsChargeSum" in histName or "cluster_neighbor_charge_sum" in histName:
        # is_target_hist = ("clusterNeighborsChargeSum" in histName) or ("cluster_neighbor_charge_sum" in histName) # この行がない
        #     # 最初のビン (X=0 または最初のビン) の内容をクリア
        #     # TAxis::GetFirst() はアンダーフロービンを除いた最初の有効なビンを返す (通常はBin 1)
        # if is_target_hist:
        #     # 除去するビンのリスト: Bin 0 (アンダーフロー) と Bin 1 (最初の有効ビン)
        #     bins_to_clear = [0, 1]

        #     for bin_index in bins_to_clear:
        #         bin_content = hist.GetBinContent(bin_index)
                
        #         if bin_content > 0:
        #             # ビンの内容をゼロに設定
        #             hist.SetBinContent(bin_index, 0)
        #             # エラーもゼロに設定
        #             hist.SetBinError(bin_index, 0) 
                    
        #             # イベント総数を減らす
        #             hist.SetEntries(hist.GetEntries() - bin_content)
                    
        #             log.info(f"getScaledHist | Cleared Bin {bin_index} content for {histName}: {bin_content:.0f} entries removed.")  
        # -----------------------------------------------------------------

        return self.scalingHistogram(hist, scaleFactor, title)

    def fitLandau(self, hist, xMin, xMax, color):
        """
        裾野（立ち上がり）の影響を排除し、ピーク位置(MPV)を正確に合わせるための修正版
        """
        func_name = f"fit_{hist.GetName()}_{ROOT.gRandom.Integer(10000)}"
        
        # ヒストグラムの基本情報を取得
        mean = hist.GetMean()
        sigma_raw = hist.GetRMS()
        max_val = hist.GetMaximum()
        max_bin = hist.GetMaximumBin()
        peak_x = hist.GetXaxis().GetBinCenter(max_bin)
        
        # --- 範囲決定のロジック変更 ---
        # 立ち上がりの裾野は「純粋なLandau」と「実データ」で形状が合わないため無視する。
        # ピーク位置(peak_x)から、左側に 0.8シグマ 〜 1.0シグマ 程度戻ったところを開始点とする。
        # これにより、立ち上がりの一番急な部分と頂点のみを使うことになる。
        
        # ここの係数(0.8)を大きくするほど、左側の裾野を無視します。
        # まだ左にズレるようなら 1.0 や 1.2 に上げてください。
        fit_min_limit = peak_x - (sigma_raw * 0.8) 
        
        # 右側は十分なデータがあるので広めにとる
        fit_max_limit = peak_x + (sigma_raw * 3.0)

        # ユーザー指定の範囲(xMin, xMax)の中に収めるガード処理
        if fit_min_limit < xMin: fit_min_limit = xMin
        if fit_max_limit > xMax: fit_max_limit = xMax
        
        # もしピーク位置推定がずれて変な範囲になった場合の安全策
        if fit_min_limit >= peak_x:
            fit_min_limit = peak_x - 100 # 強制的に少し戻す

        log.info(f"fitLandau | Range optimization: Peak={peak_x:.1f}, Range=[{fit_min_limit:.1f}, {fit_max_limit:.1f}]")

        # --- 関数定義とパラメータ設定 ---
        func = ROOT.TF1(func_name, "landau", xMin, xMax)
        
        # パラメータ初期値: [Constant, MPV, Sigma]
        # MPVの初期値はヒストグラムの最大ビン位置にする
        func.SetParameters(max_val * 5.0, peak_x, sigma_raw * 0.15)
        
        # MPV(Parameter 1)が極端にズレないように制限をかける（これ重要です）
        # ピーク位置 ± RMSの範囲内しか動かないようにする
        func.SetParLimits(1, peak_x - sigma_raw, peak_x + sigma_raw)

        # --- フィット実行 ---
        # まずは狭い範囲でフィットしてパラメータを確定させる
        func.SetRange(fit_min_limit, fit_max_limit)
        hist.Fit(func, "Q R 0 N") # Q:Quiet, R:Range, N:NoStore

        # --- 描画用調整 ---
        # フィット結果のパラメータを保持したまま、描画範囲だけ指定のxMin, xMaxに戻す
        # これにより「線」は全体に引かれるが、「形状」はピーク付近で決定されたものになる
        func.SetRange(xMin, xMax) 
        
        # スタイル設定
        func.SetLineColor(color)
        func.SetLineStyle(2) 
        func.SetLineWidth(2)

        return func
    
    def fitLangau(self, hist, xMin, xMax, color):
        """
        LandauとGaussianの畳み込み関数(Langau)でフィットを行う。
        エラー回避のため、ROOT.langaufun を直接TF1に渡す。
        """
        func_name = f"fit_langau_{hist.GetName()}_{ROOT.gRandom.Integer(10000)}"
        
        # --- パラメータ推定 ---
        max_bin = hist.GetMaximumBin()
        peak_x = hist.GetXaxis().GetBinCenter(max_bin)
        
        # Area推定
        area = hist.Integral() * hist.GetBinWidth(1)
        
        # 幅の推定
        rms = hist.GetRMS()
        landau_width_guess = rms * 0.15
        gauss_sigma_guess  = rms * 0.10

        # --- 関数定義 (ここが修正ポイント) ---
        # 文字列 "langaufun" ではなく、ROOT.langaufun を直接渡すことでリンクエラーを防ぐ
        func = ROOT.TF1(func_name, ROOT.langaufun, xMin, xMax, 4)
        
        func.SetParNames("LandauWidth", "MPV", "Area", "GausSigma")
        func.SetParameters(landau_width_guess, peak_x, area, gauss_sigma_guess)

        # --- パラメータ制限 ---
        func.SetParLimits(0, 0.1, rms)                    # Landau Width
        func.SetParLimits(1, peak_x - rms, peak_x + rms)  # MPV
        func.SetParLimits(2, area * 0.1, area * 10.0)     # Area
        func.SetParLimits(3, 1.0, rms)                    # Gaus Sigma

        # --- 描画の滑らかさ設定 ---
        func.SetNpx(500) # 描画ポイントを増やしてカクつきを防止

        # --- フィット実行 ---
        # 2回実行すると安定しやすい
        # 1回目: 範囲限定
        fit_range_min = peak_x - rms
        fit_range_max = peak_x + 3 * rms
        if fit_range_min < xMin: fit_range_min = xMin
        if fit_range_max > xMax: fit_range_max = xMax
        
        func.SetRange(fit_range_min, fit_range_max)
        hist.Fit(func, "Q R 0 N") # Pre-fit
        
        # 2回目: 全体
        func.SetRange(xMin, xMax)
        func.SetLineColor(color)
        func.SetLineStyle(2)
        func.SetLineWidth(2)
        
        # "B":境界制約あり, "M":Minuit改良
        hist.Fit(func, "Q R B M 0 N") 
        
        log.info(f"fitLangau | {hist.GetName()} MPV={func.GetParameter(1):.1f}, L_Width={func.GetParameter(0):.1f}, G_Sigma={func.GetParameter(3):.1f}")

        return func

    def calculateResidual(self, hists_info: List[Tuple[ROOT.TH1, str, str]], exp_base_key: str, data_entries: Dict) -> List[ROOT.TGraphErrors]:
        """
        指定された exp_base_key に対応するヒストグラムを基準 (Base) とし、
        他の全てのヒストグラムとの相対残差 (Comp - Base) / Base を計算し、TGraphErrorsのリストとして返す。
        """
        if not hists_info:
            return []
        
        # 基準となるExpヒストグラムを特定
        h_base_raw = None
        base_data_config = None
        
        # hists_info: [(ヒストグラム, データセット名, データキー)]
        for h_raw, name, key in hists_info:
            if key == exp_base_key:
                h_base_raw = h_raw
                base_data_config = data_entries[key]
                break
        
        if not h_base_raw:
            log.error(f"calculateResidual | Base Exp histogram not found for key: {exp_base_key}")
            return []

        # 基準Baseヒストグラムを正規化
        h_base_norm = h_base_raw.Clone(f"{h_base_raw.GetName()}_norm_base")
        h_base_norm.SetDirectory(0)
        base_integral = h_base_norm.Integral(0, h_base_norm.GetNbinsX() + 1)
        if base_integral > 0:
            h_base_norm.Scale(1.0 / base_integral)

        residuals = []
        
        # 全てのヒストグラムと比較
        for h_raw, name, key in hists_info:
            # 基準Base自体との比較はスキップ
            if key == exp_base_key:
                continue
            
            # 比較対象ヒストグラムを正規化
            h_comp_norm = h_raw.Clone(f"{h_raw.GetName()}_norm_comp")
            h_comp_norm.SetDirectory(0)
            comp_integral = h_comp_norm.Integral(0, h_comp_norm.GetNbinsX() + 1)
            if comp_integral > 0:
                h_comp_norm.Scale(1.0 / comp_integral)

            # スタイル設定を取得
            ds_comp = data_entries[key]
            color = ds_comp.color 
            marker = ds_comp.markerStyle 

            n_bins = h_base_norm.GetNbinsX()
            g_res = ROOT.TGraphErrors(n_bins)
            g_res.SetName(f"g_res_{key}")
            g_res.SetTitle(f"Residual_{name}")
            
            g_res.SetMarkerColor(color)
            g_res.SetMarkerStyle(marker)
            g_res.SetMarkerSize(0.8)
            g_res.SetLineColor(color)
            g_res.SetLineWidth(2)
            g_res.SetFillStyle(0)

            point_count = 0
            for bin_i in range(1, n_bins + 1):
                base_content = h_base_norm.GetBinContent(bin_i)
                comp_content = h_comp_norm.GetBinContent(bin_i)
                bin_center = h_base_norm.GetBinCenter(bin_i)
                
                if base_content > 0:
                    residual = (comp_content - base_content) / base_content
                    
                    # 誤差計算 (ポアソン統計を仮定した正規化後の相対誤差)
                    base_error_norm = h_base_norm.GetBinError(bin_i)
                    comp_error_norm = h_comp_norm.GetBinError(bin_i)

                    # Sim/Expの比率の誤差
                    ratio_error = (1.0 / base_content) * ROOT.TMath.Sqrt(base_error_norm**2 + comp_error_norm**2)
                    
                    g_res.SetPoint(point_count, bin_center, residual)
                    g_res.SetPointError(point_count, 0, ratio_error)
                    point_count += 1
                
            g_res.Set(point_count)
            residuals.append(g_res)
            
        return residuals

    # def plot2DComparison(self, c, l, datasets: List[DataConfig], suffix: str,
    #                      histKey: str, histLegendSuffix: str,
    #                      outName: str, plotTitle: str,
    #                      doFit: bool = False):
    #     """
    #     - doFit=True : TProfile + 2d_fitting
    #     - doFit=False: Canvas divide and plot (COLZ)
    #     """
    #     c.Clear(); l.Clear()
    #     ROOT.gStyle.SetPalette(ROOT.kBird)

    #     valid_datasets = []
    #     for ds in datasets:
    #         if ds.histNames.get(histKey):
    #             valid_datasets.append(ds)
        
    #     if not valid_datasets:
    #         log.warning(f"plot2DComparison | No datasets found for key: {histKey}")
    #         return

    #     # ---------------------------------------------------------
    #     # doFit = True
    #     # ---------------------------------------------------------
    #     if doFit:
    #         self.setCanvasMargins(c, left=0.14, right=0.06, top=0.10, bottom=0.14)
            
    #         profiles = []
    #         funcs = []
    #         max_y = 0.0

    #         for i, ds in enumerate(valid_datasets):
    #             histName = ds.histNames.get(histKey)
    #             # 2Dヒストグラム取得 (スケールは1.0推奨)
    #             h2 = self.getScaledHist(ds.fileKey, histName, 1.0, "", ds.isMerged)
                
    #             if not h2 or not h2.InheritsFrom("TH2"):
    #                 continue

    #             # 2D -> 1D Profile (X軸ビンごとのYの平均値)
    #             prof_name = f"{h2.GetName()}_prof_{ds.name}_{suffix}"
    #             prof = h2.ProfileX(prof_name)
                
    #             # スタイル設定
    #             self.setHistStyle(prof, ds.color, ds.markerStyle)
    #             prof.SetLineWidth(2)
                
    #             # フィッティング (例: pol1 = 直線近似)
    #             # ※ 必要に応じて "pol2" や ガウスなどに変更可能
    #             fit_func = ROOT.TF1(f"fit_{prof_name}", "pol1", h2.GetXaxis().GetXmin(), h2.GetXaxis().GetXmax())
    #             fit_func.SetLineColor(ds.color)
    #             fit_func.SetLineStyle(2) # 点線
                
    #             # フィット実行
    #             prof.Fit(fit_func, "Q R")
                
    #             profiles.append((prof, ds.name))
    #             funcs.append(fit_func)
                
    #             # 軸範囲調整用
    #             if prof.GetMaximum() > max_y: max_y = prof.GetMaximum()

    #         # 描画ループ (重ね書き)
    #         if not profiles: return

    #         for i, ((prof, name), func) in enumerate(zip(profiles, funcs)):
    #             opt = "PE" if i == 0 else "same PE"
                
    #             if i == 0:
    #                 # 軸タイトルなどは元の2Dヒストグラムから引き継ぐか、ここで設定
    #                 prof.SetTitle(f"{plotTitle} (Profile Fit)")
    #                 prof.GetYaxis().SetTitle(f"Average {prof.GetYaxis().GetTitle()}")
    #                 prof.GetYaxis().SetRangeUser(0, max_y * 1.5)
    #                 self.setAxisStyle(prof)
                
    #             prof.Draw(opt)
    #             func.Draw("same")
                
    #             # レジェンド (傾きなどを表示)
    #             p0 = func.GetParameter(0)
    #             p1 = func.GetParameter(1)
    #             l.AddEntry(prof, f"{name}: {p0:.2f} + {p1:.2e}x", "pl")

    #         l.Draw()
    #         self.drawTitle(f"{plotTitle} (Trend Comparison)")
    #         c.SaveAs(f"./plot/{outName}_FitCompare_{suffix}.pdf")

    #     # ---------------------------------------------------------
    #     # doFit = False
    #     # ---------------------------------------------------------
    #     else:
    #         n_plots = len(valid_datasets)
    #         # 分割数を計算 (例: 4つなら 2x2, 5つなら 3x2)
    #         n_cols = math.ceil(math.sqrt(n_plots))
    #         n_rows = math.ceil(n_plots / n_cols)
            
    #         c.Divide(n_cols, n_rows)
            
    #         for i, ds in enumerate(valid_datasets):
    #             c.cd(i + 1)
    #             ROOT.gPad.SetLeftMargin(0.12)
    #             ROOT.gPad.SetRightMargin(0.12) # カラーバー用
                
    #             histName = ds.histNames.get(histKey)
    #             h2 = self.getScaledHist(ds.fileKey, histName, 1.0, "", ds.isMerged)
                
    #             if h2 and h2.InheritsFrom("TH2"):
    #                 h2.SetTitle(f"{ds.name}")
    #                 self.setAxisStyle(h2)
                    
    #                 # 統計ボックスを消す
    #                 h2.SetStats(0)
    #                 h2.Draw("COLZ")
                    
    #                 # テキスト情報を少し入れる
    #                 t = ROOT.TLatex()
    #                 t.SetNDC()
    #                 t.SetTextSize(0.04)
    #                 t.DrawLatex(0.15, 0.85, f"{ds.name}")

    #         c.SaveAs(f"{self.output_dirs['2D']}/{outName}_Maps_{suffix}.pdf")

    def plot2DComparison(self, c, l, datasets, suffix, histKey, histLegendSuffix, outName, plotTitle, doProfile=False):
        """
        - doProfile=True: ProfileX comparison (No Fit)
        - doProfile=False: 2D Colormap
        """
        c.Clear(); l.Clear()
        ROOT.gStyle.SetPalette(ROOT.kBird)
        valid_datasets = [ds for ds in datasets if ds.histNames.get(histKey)]
        if not valid_datasets: return

        if doProfile:
            # --- Change: ProfileX Comparison without Fitting ---
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
            # 名前を FitCompare から ProfileCompare に変更
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

    # def run_plotCheck(self):
    #     """メイン実行関数"""
    #     log.info("run_plotCheck | Start")

    #     data_entries = load_data_config(self.config_file)

    #     if not data_entries:
    #         log.critical("No data entries loaded. Please check config.json.")
    #         return

    #     # --- データセット定義 ---
    #     charge_rebin_exp_defaults = {
    #         "clusterCharge": 20,
    #         "seedCharge": 20,
    #         "neighborChargeSum": 10,
    #         "seedChargeBySize": 20, 
    #         "clusterSize": 1,
    #         "driftTime": 1,
    #         "default": 1
    #     }
    #     charge_rebin_sim_defaults = {
    #         "clusterCharge": 4,
    #         "seedCharge": 4,
    #         "neighborChargeSum": 1,
    #         "seedChargeBySize": 4, 
    #         "clusterSize": 1,
    #         "driftTime": 1,
    #         "default": 1
    #     }

    #     for key,config in data_entries.items():
    #         if "sim" in key.lower() or "sim" in config.name.lower():
    #             defaults = charge_rebin_sim_defaults.copy()
    #         else:
    #             defaults = charge_rebin_exp_defaults.copy()
            
    #         defaults.update(config.rebinValues)
    #         config.rebinValues = defaults
            
    #         log.info(f"Applied rebin settings for {config.name}: {config.rebinValues}")

    #     comparisons = [
    #         # # 1. GAPの Exp vs Sim比較
    #         # #(["exp_gap", "exp_gap_nt600", "sim_gap_masetti_et2", "sim_gap_masetti_et2_nt144", "sim_gap_masetti_et2_nt48_n90e"], "_gap_default"),
    #         # (["exp_gap_nt600", "sim_gap_masetti_et2_nt144_n90e"], "GAP_nt600"),
    #         # #(["exp_gap_nt600", "sim_gap_masetti_et2_nt48_n90e_srh_dortmund", "sim_gap_canali_et2_nt48_n90e_srh_dortmund", "sim_gap_jacoboni_et2_nt48_n90e_srh_dortmund"], "GAP_CMModel"),
    #         # (["exp_gap", "sim_gap_masetti_et2_nt48_n90e_srh_dortmund", "sim_gap_canali_et2_nt48_n90e_srh_dortmund", "sim_gap_jacoboni_et2_nt48_n90e_srh_dortmund"], "GAP_CMModel"),
    #         # #(["exp_gap", "sim_gap_jacoboni_et2_nt48_n90e_srh_dortmund"], "_gap_nt48_CMmodel"),
    #         # (["exp_gap", "sim_gap_masetti_et2_nt48_n90e_srh_withoutTP", "sim_gap_masetti_et2_nt48_n90e_auger_withoutTP", "sim_gap_masetti_et2_nt48_n90e_srhauger_withoutTP"], "GAP_RCModel"),
    #         # (["exp_gap", "sim_gap_masetti_et2_nt48_n90e_srh_mandic", "sim_gap_masetti_et2_nt48_n90e_srh_ljubljana", "sim_gap_masetti_et2_nt48_n90e_srh_cmstracker", "sim_gap_masetti_et2_nt48_n90e_srh_dortmund"], "GAP_TPModel"),
    #         # #(["exp_gap", "sim_gap_masetti_et2"], "_gap_default"),
    #         # #(["exp_gap_nt600", "sim_gap_masetti_et2_nt120", "sim_gap_masetti_et2_nt144"], "_gap_default"),
    #         (["exp_p15_gap_10v_st1000_nt600",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_jacoboni_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_canali_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburg_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_masetti_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_arora_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_Ecanali_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_quay_none_none"
    #           ], "p15_GAP_10V_st240_nt144_CMModel"),
    #           (["exp_p15_gap_10v_st1000_nt600",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_jacoboni_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_canali_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburg_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_quay_none_none",
    #           ], "p15_GAP_efield"),
    #           (["exp_p15_gap_10v_st1000_nt600",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_masetti_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_arora_none_none",
    #           "sim_p15_gap_10v_n90e_303K_st240_nt144_Ecanali_none_none"
    #           ], "p15_GAP_doping"),


    #           (["exp_p22p5_gap_10v_st1000_nt600",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_jacoboni_none_none",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_canali_none_none",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_hamburg_none_none",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_quay_none_none",
    #           ], "p22p5_GAP_efield"),
    #           (["exp_p22p5_gap_10v_st1000_nt600",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_masetti_none_none",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_arora_none_none",
    #           "sim_p22p5_gap_10v_n90e_303K_st240_nt144_Ecanali_none_none"
    #           ], "p22p5_GAP_doping"),
           
            

    #         # # 2. STDの Exp vs Sim比較
    #         # #(["exp_std", "exp_std_nt600", "sim_std_masetti_et2", "sim_std_masetti_et2_nt144", "sim_std_masetti_et2_nt48_n90e"], "_std_default"),
    #         # (["exp_std_nt600","sim_std_masetti_et2_nt144_n90e"], "STD_nt600"),
    #         # (["exp_std", "sim_std_masetti_et2_nt48_n90e_srh_dortmund", "sim_std_canali_et2_nt48_n90e_srh_dortmund", "sim_std_jacoboni_et2_nt48_n90e_srh_dortmund"], "STD_CMModel"),
    #         # #(["exp_std", "sim_std_jacoboni_et2_nt48_n90e_srh_dortmund"], "_std_nt48_CMmodel"),
    #         # (["exp_std", "sim_std_masetti_et2_nt48_n90e_srh_withoutTP", "sim_std_masetti_et2_nt48_n90e_auger_withoutTP", "sim_std_masetti_et2_nt48_n90e_srhauger_withoutTP"], "STD_RCModel"),
    #         # (["exp_std", "sim_std_masetti_et2_nt48_n90e_srh_mandic", "sim_std_masetti_et2_nt48_n90e_srh_ljubljana", "sim_std_masetti_et2_nt48_n90e_srh_cmstracker", "sim_std_masetti_et2_nt48_n90e_srh_dortmund"], "STD_TPModel"),
    #         # #(["exp_std", "sim_std_masetti_et2"], "_std_default"),
    #         # #(["exp_std_nt600", "sim_std_masetti_et2_nt120", "sim_std_masetti_et2_nt144"], "_std_default"),
    #         (["exp_p15_std_10v_st1000_nt600",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_jacoboni_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_canali_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_hamburg_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_masetti_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_arora_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_Ecanali_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_quay_none_none"
    #           ], "p15_STD_10V_st240_nt144_CMModel"),
    #           (["exp_p15_std_10v_st1000_nt600",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_jacoboni_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_canali_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_hamburg_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_quay_none_none",
    #           ], "p15_STD_efield"),
    #           (["exp_p15_std_10v_st1000_nt600",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_masetti_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_arora_none_none",
    #           "sim_p15_std_10v_n90e_303K_st240_nt144_Ecanali_none_none",
    #           ], "p15_STD_doping"),


    #           (["exp_p22p5_std_10v_st1000_nt600",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_jacoboni_none_none",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_canali_none_none",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_hamburg_none_none",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_quay_none_none",
    #           ], "p22p5_STD_efield"),

    #           (["exp_p22p5_std_10v_st1000_nt600",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_masetti_none_none",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_arora_none_none",
    #           "sim_p22p5_std_10v_n90e_303K_st240_nt144_Ecanali_none_none",
    #           ], "p22p5_STD_doping")
    #     ]

    #     canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    #     legend = ROOT.TLegend(0.35, 0.55, 0.9, 0.8)
    #     legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.03)

    #     for keys, suffix in comparisons:
    #         datasets = [data_entries[k] for k in keys if k in data_entries]
    #         if not datasets: continue

    #         log.info(f"Processing comparison set: {keys}")
            
    #         plot_title = " vs ".join([d.name for d in datasets])
    #         if len(plot_title) > 50: plot_title = f"{suffix}"

    #         self.plotAllCharge(canvas, legend, datasets, suffix, plot_title)
    #         self.plotClusterSize(canvas, legend, datasets, suffix, plot_title)

    #         self.plotComparison(canvas, legend, datasets, suffix,
    #             histKey="seedCharge", histLegendSuffix="seed",
    #             outName="clusterSeedCharge_comp", plotTitle=plot_title,
    #             xMin=0, xMax=4000,
    #             doLandauFit=False) 
            
    #         self.plotComparison(canvas, legend, datasets, suffix,
    #             histKey="clusterCharge", histLegendSuffix="cluster",
    #             outName="clusterCharge_comp", plotTitle=plot_title,
    #             xMin=0, xMax=4000,
    #             doLandauFit=True) # ClusterChargeはFitする

    #         self.plotComparison(canvas, legend, datasets, suffix,
    #             histKey="neighborChargeSum", histLegendSuffix="neighbor sum",
    #             outName="clusterNeighborChargeSum_comp", plotTitle=plot_title,
    #             xMin=20, xMax=2000,
    #             doLandauFit=False)
            
    #         self.plotComparison(canvas, legend, datasets, suffix,
    #             histKey="driftTime", histLegendSuffix="drift",
    #             outName="driftTime_comp", plotTitle=plot_title,
    #             xMin=0, xMax=0.2,
    #             doLandauFit=False)

    #         self.plotSeedChargeByCS(canvas, legend, datasets, suffix, plot_title)

    #         self.plot2DComparison(canvas, legend, datasets, suffix,
    #                               histKey="clusterSizeVsCharge", # JSONに追加したキー
    #                               histLegendSuffix="Size vs Charge",
    #                               outName="2D_SizeVsCharge",
    #                               plotTitle=plot_title,
    #                               doFit=True)
            
    #         c_div = ROOT.TCanvas("c_div", "c_div", 1200, 800)
    #         self.plot2DComparison(c_div, legend, datasets, suffix,
    #                               histKey="clusterSizeVsCharge",
    #                               histLegendSuffix="Size vs Charge",
    #                               outName="2D_SizeVsCharge",
    #                               plotTitle=plot_title,
    #                               doFit=False)
        
    #         canvas.cd()

    #         # === 追加: Path Analysis (2D map + 1D graph per dataset) ===
    #         #  "clusterSizeVsCharge" キーなど、2Dマップが定義されているものを指定
    #         #  JSON側で 2Dヒストグラムのキー (例: "efficiencyMap" や "residualMap") があると想定
            
    #         target_2d_key = "clusterSizeVsCharge" # ※適切な2Dヒストグラムのキーに変更してください
            
    #         # 1. 個別データセットごとのプロット
    #         for ds in datasets:
    #             # このデータセットに2Dマップの設定があるか確認
    #             if target_2d_key in ds.histNames:
    #                 self.plotPathAnalysis(ds, target_2d_key, suffix)

    #         # 2. 結合プロット (比較)
    #         self.plotPathAnalysisCombined(canvas, legend, datasets, target_2d_key, suffix, plot_title)
        
    #     self.cleanupResources()
    #     log.info("run_plotCheck | Finished")

    def run_plotCheck(self):
        """メイン実行関数"""
        log.info("run_plotCheck | Start")

        data_entries = load_data_config(self.config_file)

        if not data_entries:
            log.critical("No data entries loaded. Please check config.json.")
            return

        # --- データセット定義 ---
        charge_rebin_exp_defaults = {
            "clusterCharge": 20,
            "seedCharge": 20,
            "neighborChargeSum": 10,
            "seedChargeBySize": 20, 
            "clusterSize": 1,
            "driftTime": 1,
            "default": 1
        }
        charge_rebin_sim_defaults = {
            "clusterCharge": 4,
            "seedCharge": 4,
            "neighborChargeSum": 1,
            "seedChargeBySize": 4, 
            "clusterSize": 1,
            "driftTime": 1,
            "default": 1
        }

        for key, config in data_entries.items():
            if "sim" in key.lower() or "sim" in config.name.lower():
                defaults = charge_rebin_sim_defaults.copy()
            else:
                defaults = charge_rebin_exp_defaults.copy()
            
            defaults.update(config.rebinValues)
            config.rebinValues = defaults
            
            log.info(f"Applied rebin settings for {config.name}: {config.rebinValues}")

        comparisons = [
            # # 1. GAPの Exp vs Sim比較
            # #(["exp_gap", "exp_gap_nt600", "sim_gap_masetti_et2", "sim_gap_masetti_et2_nt144", "sim_gap_masetti_et2_nt48_n90e"], "_gap_default"),
            # (["exp_gap_nt600", "sim_gap_masetti_et2_nt144_n90e"], "GAP_nt600"),
            # #(["exp_gap_nt600", "sim_gap_masetti_et2_nt48_n90e_srh_dortmund", "sim_gap_canali_et2_nt48_n90e_srh_dortmund", "sim_gap_jacoboni_et2_nt48_n90e_srh_dortmund"], "GAP_CMModel"),
            # (["exp_gap", "sim_gap_masetti_et2_nt48_n90e_srh_dortmund", "sim_gap_canali_et2_nt48_n90e_srh_dortmund", "sim_gap_jacoboni_et2_nt48_n90e_srh_dortmund"], "GAP_CMModel"),
            # #(["exp_gap", "sim_gap_jacoboni_et2_nt48_n90e_srh_dortmund"], "_gap_nt48_CMmodel"),
            # (["exp_gap", "sim_gap_masetti_et2_nt48_n90e_srh_withoutTP", "sim_gap_masetti_et2_nt48_n90e_auger_withoutTP", "sim_gap_masetti_et2_nt48_n90e_srhauger_withoutTP"], "GAP_RCModel"),
            # (["exp_gap", "sim_gap_masetti_et2_nt48_n90e_srh_mandic", "sim_gap_masetti_et2_nt48_n90e_srh_ljubljana", "sim_gap_masetti_et2_nt48_n90e_srh_cmstracker", "sim_gap_masetti_et2_nt48_n90e_srh_dortmund"], "GAP_TPModel"),
            # #(["exp_gap", "sim_gap_masetti_et2"], "_gap_default"),
            # #(["exp_gap_nt600", "sim_gap_masetti_et2_nt120", "sim_gap_masetti_et2_nt144"], "_gap_default"),
            (["exp_p15_gap_10v_st1000_nt600",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_jacoboni_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_canali_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburg_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_masetti_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_arora_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_Ecanali_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_quay_none_none"
              ], "p15_GAP_10V_st240_nt144_CMModel"),
              (["exp_p15_gap_10v_st1000_nt600",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_jacoboni_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_canali_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburg_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_quay_none_none",
              ], "p15_GAP_efield"),
              (["exp_p15_gap_10v_st1000_nt600",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_masetti_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_arora_none_none",
              "sim_p15_gap_10v_n90e_303K_st240_nt144_Ecanali_none_none"
              ], "p15_GAP_doping"),


              (["exp_p22p5_gap_10v_st1000_nt600",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_jacoboni_none_none",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_canali_none_none",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_hamburg_none_none",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_quay_none_none",
              ], "p22p5_GAP_efield"),
              (["exp_p22p5_gap_10v_st1000_nt600",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_masetti_none_none",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_arora_none_none",
              "sim_p22p5_gap_10v_n90e_303K_st240_nt144_Ecanali_none_none"
              ], "p22p5_GAP_doping"),
           
            

            # # 2. STDの Exp vs Sim比較
            # #(["exp_std", "exp_std_nt600", "sim_std_masetti_et2", "sim_std_masetti_et2_nt144", "sim_std_masetti_et2_nt48_n90e"], "_std_default"),
            # (["exp_std_nt600","sim_std_masetti_et2_nt144_n90e"], "STD_nt600"),
            # (["exp_std", "sim_std_masetti_et2_nt48_n90e_srh_dortmund", "sim_std_canali_et2_nt48_n90e_srh_dortmund", "sim_std_jacoboni_et2_nt48_n90e_srh_dortmund"], "STD_CMModel"),
            # #(["exp_std", "sim_std_jacoboni_et2_nt48_n90e_srh_dortmund"], "_std_nt48_CMmodel"),
            # (["exp_std", "sim_std_masetti_et2_nt48_n90e_srh_withoutTP", "sim_std_masetti_et2_nt48_n90e_auger_withoutTP", "sim_std_masetti_et2_nt48_n90e_srhauger_withoutTP"], "STD_RCModel"),
            # (["exp_std", "sim_std_masetti_et2_nt48_n90e_srh_mandic", "sim_std_masetti_et2_nt48_n90e_srh_ljubljana", "sim_std_masetti_et2_nt48_n90e_srh_cmstracker", "sim_std_masetti_et2_nt48_n90e_srh_dortmund"], "STD_TPModel"),
            # #(["exp_std", "sim_std_masetti_et2"], "_std_default"),
            # #(["exp_std_nt600", "sim_std_masetti_et2_nt120", "sim_std_masetti_et2_nt144"], "_std_default"),
            (["exp_p15_std_10v_st1000_nt600",
              "sim_p15_std_10v_n90e_303K_st240_nt144_jacoboni_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_canali_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_hamburg_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_masetti_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_arora_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_Ecanali_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_quay_none_none"
              ], "p15_STD_10V_st240_nt144_CMModel"),
              (["exp_p15_std_10v_st1000_nt600",
              "sim_p15_std_10v_n90e_303K_st240_nt144_jacoboni_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_canali_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_hamburg_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_quay_none_none",
              ], "p15_STD_efield"),
              (["exp_p15_std_10v_st1000_nt600",
              "sim_p15_std_10v_n90e_303K_st240_nt144_masetti_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_arora_none_none",
              "sim_p15_std_10v_n90e_303K_st240_nt144_Ecanali_none_none",
              ], "p15_STD_doping"),


              (["exp_p22p5_std_10v_st1000_nt600",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_jacoboni_none_none",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_canali_none_none",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_hamburg_none_none",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_hamburgHF_none_none",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_quay_none_none",
              ], "p22p5_STD_efield"),

              (["exp_p22p5_std_10v_st1000_nt600",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_masetti_none_none",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_arora_none_none",
              "sim_p22p5_std_10v_n90e_303K_st240_nt144_Ecanali_none_none",
              ], "p22p5_STD_doping")
        ]

        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        legend = ROOT.TLegend(0.35, 0.55, 0.9, 0.8)
        legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.03)

        for keys, suffix in comparisons:
            datasets = [data_entries[k] for k in keys if k in data_entries]
            if not datasets: continue
            log.info(f"Processing set: {keys}")
            plot_title = " vs ".join([d.name for d in datasets])
            if len(plot_title) > 50: plot_title = f"{suffix}"

            self.plotAllCharge(canvas, legend, datasets, suffix, plot_title)
            self.plotClusterSize(canvas, legend, datasets, suffix, plot_title)
            self.plotComparison(canvas, legend, datasets, suffix, "seedCharge", "seed", "clusterSeedCharge_comp", plot_title, 0, 4000, False) 
            self.plotComparison(canvas, legend, datasets, suffix, "clusterCharge", "cluster", "clusterCharge_comp", plot_title, 0, 4000, True)
            self.plotComparison(canvas, legend, datasets, suffix, "neighborChargeSum", "neighbor sum", "clusterNeighborChargeSum_comp", plot_title, 20, 2000, False)
            self.plotComparison(canvas, legend, datasets, suffix, "driftTime", "drift", "driftTime_comp", plot_title, 0, 0.2, False)
            self.plotSeedChargeByCS(canvas, legend, datasets, suffix, plot_title)

            # --- Change: Updated 2D plot config with PathAnalysis flag ---
            # (histKey, LegendSuffix, OutName, doProfile, doPathAnalysis)
            map_plots_config = [
                # Spatial Maps -> doPathAnalysis=True
                ("map_seedCharge", "Seed Charge Map", "2D_Map_SeedCharge", False, True),
                ("map_clusterSize", "Cluster Size Map", "2D_Map_ClusterSize", False, True),
                ("map_neighborChargeSum", "Neighbor Sum Map", "2D_Map_NeighborSum", False, True),
                
                # Correlation Maps -> doProfile=True (ProfileX), doPathAnalysis=False
                ("map_clusterVSCSize", "Cluster Charge vs Size", "2D_Corr_ClusterVSSize", True, False),
                ("map_seedVSCSize", "Seed Charge vs Size", "2D_Corr_SeedVSSize", True, False),
                ("map_seedVSneighborChargeSum", "Seed vs Neighbor Sum", "2D_Corr_SeedVSNeighbor", False, False),
            ]

            for hist_key, leg_suffix, out_name, do_profile, do_path_analysis in map_plots_config:
                if any(ds.histNames.get(hist_key) for ds in datasets):
                    
                    # 1. Comparison Plot (Profile or Path)
                    if do_profile:
                        # Correlation Plot -> ProfileX Comparison (without Fit)
                        self.plot2DComparison(canvas, legend, datasets, suffix, hist_key, leg_suffix, out_name, plot_title, True)
                    elif do_path_analysis:
                        # Spatial Map -> Path Analysis Comparison
                        self.plotPathAnalysisCombined(canvas, legend, datasets, hist_key, suffix, plot_title)

                    # 2. Individual Map Plots (Color Map)
                    c_div = ROOT.TCanvas("c_div", "c_div", 1200, 800)
                    self.plot2DComparison(c_div, legend, datasets, suffix, hist_key, leg_suffix, out_name, plot_title, False)
                    c_div.Close()
                    
                    # 3. Individual Path Analysis (Graph + Map) - Only for spatial maps
                    if do_path_analysis:
                        for ds in datasets:
                            if ds.histNames.get(hist_key):
                                self.plotPathAnalysis(ds, hist_key, suffix)

            canvas.cd()
        
        self.cleanupResources()
        log.info("run_plotCheck | Finished")

    def setCanvasMargins(self, canvas, left=0.14, right=0.06, top=0.10, bottom=0.14):
        """
        キャンバスの余白(マージン)を設定してプロット位置を調整する
        0.0〜1.0 の比率で指定 (例: 0.1 = 10%)
        """
        if not canvas: return
        canvas.SetLeftMargin(left)    # 左余白 (Y軸タイトル用)
        canvas.SetRightMargin(right)  # 右余白
        canvas.SetTopMargin(top)      # 上余白 (タイトル用)
        canvas.SetBottomMargin(bottom)# 下余白 (X軸タイトル用)

    # --- プロット関数群 (変更なし部分は省略せず記載) ---

    def plotAllCharge(self, c, l, datasets: List[DataConfig], suffix: str, plotTitle: str):
        c.Clear(); l.Clear()
        self.setCanvasMargins(c, left=0.14, right=0.06, top=0.10, bottom=0.14)

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

                self.setAxisStyle(h_cl)
                self.setAxisStyle(h_sd)

                # h_cl.Scale(1.0 / h_cl.GetEntries())
                # h_sd.Scale(1.0 / h_sd.GetEntries())
                # if h_cl.GetMaximum() > 0: h_cl.Scale(1.0 / h_cl.GetMaximum())
                # if h_sd.GetMaximum() > 0: h_sd.Scale(1.0 / h_sd.GetMaximum())

                # if h_cl.Integral(0, h_cl.GetNbinsX() + 1) > 0: 
                #     h_cl.Scale(1.0 / h_cl.Integral(0, h_cl.GetNbinsX() + 1))
                # if h_sd.Integral(0, h_sd.GetNbinsX() + 1) > 0: 
                #     h_sd.Scale(1.0 / h_sd.Integral(0, h_sd.GetNbinsX() + 1))

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
        #c.SaveAs(f"./plot/AllclusterCharge_{suffix}.pdf")
        c.SaveAs(f"{self.output_dirs['1D']}/AllclusterCharge_{suffix}.pdf")

    def plotClusterSize(self, c, l, datasets: List[DataConfig], suffix: str, plotTitle: str):
        c.Clear(); l.Clear()
        self.setCanvasMargins(c, left=0.14, right=0.06, top=0.10, bottom=0.14)

        hists = []
        max_y = 0.0
        for ds in datasets:
            h = self.getScaledHist(ds.fileKey, ds.histNames.get("clusterSize"), 1.0, "", ds.isMerged) 
            if h:
                self.setHistStyle(h, ds.color, ds.markerStyle, 0.2)
                # if h.GetEntries() > 0: h.Scale(1.0 / h.GetEntries())
                # h.Scale(1.0 / h.GetMaximum())
                if h.Integral(0, h.GetNbinsX() + 1) > 0: 
                    h.Scale(1.0 / h.Integral(0, h.GetNbinsX() + 1))

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
            self.setAxisStyle(h)
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
        self.setCanvasMargins(c, left=0.14, right=0.06, top=0.10, bottom=0.14)

        hists = []
        funcs = [] 
        max_y = 0.0

        # Title
        title = ";charge [e];counts"
        if histKey == "driftTime":
            title = ";drift time [ns];counts"

        for ds in datasets:
            histName = ds.histNames.get(histKey)
            if not histName: continue

            h = self.getScaledHist(ds.fileKey, histName, ds.scale, "", ds.isMerged)
            if h:
                self.setHistStyle(h, ds.color, ds.markerStyle)
                rebin_val = ds.rebinValues.get(histKey, ds.rebinValues.get("default", 1))
                if rebin_val > 1 and ds.doRebin:
                    h.Rebin(rebin_val)
                if h.GetMaximum() > 0:
                    # h.Scale(1.0 / h.GetEntries())
                    h.Scale(1.0 / h.GetMaximum())
                # if h.Integral(0, h.GetNbinsX() + 1) > 0: 
                #     h.Scale(1.0 / h.Integral(0, h.GetNbinsX() + 1))

                max_y = max(max_y, h.GetMaximum())

                f = None
                if doLandauFit:
                    # ここでフィットを実行
                    #f = self.fitLandau(h, xMin, xMax, ds.color)
                    f = self.fitLangau(h, xMin, xMax, ds.color)
                    funcs.append(f)
                else:
                    funcs.append(None)
                hists.append((h, ds.name))

        if not hists: return

        for i, (h, name) in enumerate(hists):
            #opt = "PE" if i == 0 else "same PE"
            #opt = "HIST" if i == 0 else "HIST SAME"
            # h.SetMarkerSize(0)

            # opt_fill = "E3" if i == 0 else "E3 SAME"
            # opt_line = "L" if i == 0 else "L SAME"

            # # if i == 0:
            # #     #h.SetTitle(";charge [e];counts")
            # #     h.SetTitle(title)
            # #     h.GetXaxis().SetRangeUser(xMin, xMax)
            # #     h.GetYaxis().SetRangeUser(0, max_y * 1.2)
            # # self.setAxisStyle(h)
            # # #h.Draw(opt)
            # # h.Draw(opt_fill)
            # # h.Draw(opt_line + " SAME")
            
            # if i == 0:
            #     h.SetTitle(title)
            #     h.GetXaxis().SetRangeUser(xMin, xMax)
            #     h.GetYaxis().SetRangeUser(0, max_y * 1.3)
            #     self.setAxisStyle(h)
            #     h.Draw(opt_fill) # 塗りつぶし（下層）
            # else:
            #     h.Draw(opt_fill)
            
            # h.Draw(opt_line) # 線（上層）

            opt_band = "E3" if i == 0 else "E3 SAME"
            # L: データ点を結ぶ線を描画 (C に変えると滑らかな曲線になります)
            opt_line = "L SAME"
        
            if i == 0:
                h.SetTitle(title)
                h.GetXaxis().SetRangeUser(xMin, xMax)
                h.GetYaxis().SetRangeUser(0, max_y * 1.3)
                self.setAxisStyle(h)
                h.Draw(opt_band)
            else:
                h.Draw(opt_band)
            

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
        c.SaveAs(f"{self.output_dirs['1D']}/{outName}_{suffix}.pdf")

    # def plotSeedChargeByCS(self, c, l, datasets: List[DataConfig], suffix: str, plotTitleBase: str):
    #     valid_datasets = [ds for ds in datasets if ds.histNames.get("seedChargeBySize")]
    #     if not valid_datasets: return

    #     for cs in range(1, 7):
    #         c.Clear(); l.Clear()
    #         hists = []
    #         max_y = 0.0
    #         for ds in valid_datasets:
    #             baseName = ds.histNames.get("seedChargeBySize")
    #             targetName = f"{baseName}{cs}"
    #             h = self.getScaledHist(ds.fileKey, targetName, ds.scale, "", ds.isMerged)
    #             if h:
    #                 self.setHistStyle(h, ds.color, ds.markerStyle)

    #                 rebin_val = ds.rebinValues.get("seedChargeBySize", 1)
    #                 if rebin_val > 1 and ds.doRebin:
    #                     h.Rebin(rebin_val)
    #                 if h.GetMaximum() > 0: h.Scale(1.0 / h.GetEntries())
    #                 # if h.Integral(0, h.GetNbinsX() + 1) > 0: 
    #                 #     h.Scale(1.0 / h.Integral(0, h.GetNbinsX() + 1))
    #                 max_y = max(max_y, h.GetMaximum())
    #                 hists.append((h, ds.name))

    #         if not hists: continue
    #         for i, (h, name) in enumerate(hists):
    #             opt = "PE" if i == 0 else "same PE"
    #             if i == 0:
    #                 h.SetTitle(";charge [e];counts")
    #                 h.GetXaxis().SetRangeUser(0, 4000)
    #                 h.GetYaxis().SetRangeUser(0, max_y * 1.2)
    #             h.Draw(opt)
    #             l.AddEntry(h, f"{name}, seed", "pe")
    #         l.Draw()
    #         self.drawTitle(f"{plotTitleBase} cs{cs}")
    #         c.SaveAs(f"./plot/seed_charge_cs/clusterSeedCharge_cs{cs}_{suffix}.pdf")

    def plotSeedChargeByCS(self, c, l, datasets: List[DataConfig], suffix: str, plotTitleBase: str):
        # seedChargeBySize が設定されているデータセットのみ対象
        valid_datasets = [ds for ds in datasets if ds.histNames.get("seedChargeBySize")]
        if not valid_datasets: return

        for cs in range(1, 7):
            c.Clear(); l.Clear()
            self.setCanvasMargins(c, left=0.14, right=0.06, top=0.10, bottom=0.14)

            hists = []
            max_y = 0.0
            
            for ds in valid_datasets:
                baseName = ds.histNames.get("seedChargeBySize")
                
                # --- ★ここから修正: パス生成ロジックの分岐 ---
                if "Sim" in ds.name:
                    # シミュレーションの場合: フォルダ階層を含めたパスを作成
                    # 例: PerClusterSize/clsize_1/seed_charge_size_1
                    # ※ baseName が "seed_charge_size_" であることを前提としています
                    # もし baseName を使わず完全固定にするなら:
                    # targetName = f"PerClusterSize/clsize_{cs}/seed_charge_size_{cs}"
                    targetName = f"PerClusterSize/clsize_{cs}/seed_charge_size_{cs}"
                else:
                    # 実験データの場合: 従来どおり末尾に数字を結合
                    # 例: .../clusterSeedCharge_size1
                    targetName = f"{baseName}{cs}"
                # ---------------------------------------------

                # 読み込み確認用ログ（必要なければコメントアウト）
                # log.info(f"Looking for hist: {targetName} in {ds.name}")

                h = self.getScaledHist(ds.fileKey, targetName, ds.scale, "", ds.isMerged)
                
                if h:
                    self.setHistStyle(h, ds.color, ds.markerStyle)

                    rebin_val = ds.rebinValues.get("seedChargeBySize", 1)
                    if rebin_val > 1 and ds.doRebin:
                        h.Rebin(rebin_val)
                    
                    #if h.GetMaximum() > 0: h.Scale(1.0 / h.GetEntries())
                    if h.GetMaximum() > 0: h.Scale(1.0 / h.GetMaximum())
                    
                    max_y = max(max_y, h.GetMaximum())
                    hists.append((h, ds.name))
                else:
                    # 見つからなかった場合の警告
                    log.warning(f"Hist not found: {targetName} in {ds.name}")

            if not hists: continue
            
            for i, (h, name) in enumerate(hists):
                opt = "PE" if i == 0 else "same PE"
                if i == 0:
                    h.SetTitle(f"Seed Charge (Size {cs});charge [e];counts")
                    h.GetXaxis().SetRangeUser(0, 4000)
                    h.GetYaxis().SetRangeUser(0, max_y * 1.2)
                h.Draw(opt)
                l.AddEntry(h, f"{name}, seed", "pe")
            
            l.Draw()
            self.drawTitle(f"{plotTitleBase} size {cs}")
            c.SaveAs(f"./plot/seed_charge_cs/clusterSeedCharge_cs{cs}_{suffix}.pdf")

    def drawTitle(self, text):
        t = ROOT.TLatex()
        t.SetTextAlign(12); t.SetTextSize(0.05)
        t.DrawLatexNDC(0.5, 0.85, text)

    # ------------------------------------------------------------------------
    #  In-Pixel Path Analysis Helper
    # ------------------------------------------------------------------------
    def extractDataAlongPath(self, hist, x1, y1, x2, y2, dist_offset, n_steps=100):
        """
        2Dヒストグラム上の (x1, y1) から (x2, y2) への線分上のデータを抽出して
        TGraphErrors を作成する。
        """
        if not hist: return None
        
        x_vals = []
        y_vals = []
        x_errs = []
        y_errs = []
        
        # 始点と終点の距離
        segment_len = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        for i in range(n_steps + 1):
            t = i / float(n_steps)
            curr_x = x1 + t * (x2 - x1)
            curr_y = y1 + t * (y2 - y1)
            
            # グローバルビン番号を取得
            bin_idx = hist.FindBin(curr_x, curr_y)
            
            # 同じビンを重複して取得しないための簡易チェック (隣接ステップと比較)
            # ※厳密には座標ごとの補間はせず、ビンの値をそのまま読む
            if i > 0:
                prev_x = x1 + (i - 1) / float(n_steps) * (x2 - x1)
                prev_y = y1 + (i - 1) / float(n_steps) * (y2 - y1)
                prev_bin = hist.FindBin(prev_x, prev_y)
                if bin_idx == prev_bin:
                    continue

            content = hist.GetBinContent(bin_idx)
            error = hist.GetBinError(bin_idx)
            
            # 経路に沿った距離
            curr_dist = dist_offset + t * segment_len
            
            x_vals.append(curr_dist)
            y_vals.append(content)
            x_errs.append(0.0)
            y_errs.append(error)
            
        if not x_vals:
            return None
            
        import array
        graph = ROOT.TGraphErrors(len(x_vals), 
                                  array.array('d', x_vals), 
                                  array.array('d', y_vals),
                                  array.array('d', x_errs), 
                                  array.array('d', y_errs))
        return graph

    # ------------------------------------------------------------------------
    #  Individual Path Analysis (2D Map + 1D Path Plot)
    # ------------------------------------------------------------------------
    def plotPathAnalysis(self, ds, histKey, suffix):
        # ... (前半部分は変更なし) ...
        histName = ds.histNames.get(histKey)
        if not histName: return None
        h2 = self.getScaledHist(ds.fileKey, histName, 1.0, "", ds.isMerged)
        if not h2: return None
        
        c_path = ROOT.TCanvas(f"c_path_{ds.name}", "Path Analysis", 1200, 600)
        c_path.Divide(2, 1)
        
        # --- Pad 1: Map ---
        pad1 = c_path.cd(1); pad1.SetRightMargin(0.15); pad1.SetLeftMargin(0.12)
        h2.SetStats(0); h2.SetTitle(f"{ds.name} 2D Map"); self.setAxisStyle(h2)
        h2.Draw("COLZ")
        
        x_max_axis = h2.GetXaxis().GetXmax(); y_max_axis = h2.GetYaxis().GetXmax()
        inset_x = x_max_axis * 0.1; inset_y = y_max_axis * 0.1
        
        pA = (0.0, 0.0)
        pB = (0.0, y_max_axis - inset_y)
        pC = (x_max_axis - inset_x, y_max_axis - inset_y)
        pD = (x_max_axis - inset_x, 0.0)
        
        # ... (矢印・マーカー描画部分は、TArrowなどをリスト管理しているためOK) ...
        # ... (省略) ...
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

        # --- Pad 2: Graph ---
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
        
        # ★ オブジェクト保持用リスト
        kept_bands = []

        for i, g in enumerate(graphs):
            if g:
                g.SetLineColor(colors[i]); g.SetLineWidth(2); g.SetMarkerColor(colors[i]); g.SetMarkerStyle(20); g.SetMarkerSize(0.5)
                g_band = g.Clone(); g_band.SetFillColorAlpha(colors[i], 0.3); g_band.SetFillStyle(1001)
                
                mg.Add(g_band, "3"); mg.Add(g, "L"); legend.AddEntry(g, labels[i], "l")
                
                # ★ リストに追加
                kept_bands.append(g_band)
        
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

    # ------------------------------------------------------------------------
    #  Combined Path Analysis (Compare multiple datasets on 1D Path)
    # ------------------------------------------------------------------------
    def plotPathAnalysisCombined(self, c, l, datasets, histKey, suffix, plotTitle):
        c.Clear(); l.Clear(); self.setCanvasMargins(c)
        mg_combined = ROOT.TMultiGraph()
        mg_combined.SetTitle(f"{plotTitle};Distance along path [um];Value")
        has_data = False
        
        # ★ オブジェクトが消えないように保持するリスト
        kept_graphs = []
        kept_bands = []

        for ds in datasets:
            histName = ds.histNames.get(histKey)
            if not histName: continue
            h2 = self.getScaledHist(ds.fileKey, histName, 1.0, "", ds.isMerged)
            if not h2: continue
            
            # --- Geometry starting from 0,0 ---
            x_max_axis = h2.GetXaxis().GetXmax(); y_max_axis = h2.GetYaxis().GetXmax()
            inset_x = x_max_axis * 0.1; inset_y = y_max_axis * 0.1
            pA = (0.0, 0.0)
            pB = (0.0, y_max_axis - inset_y)
            pC = (x_max_axis - inset_x, y_max_axis - inset_y)
            pD = (x_max_axis - inset_x, 0.0)
            
            dist_AB = math.sqrt((pB[0]-pA[0])**2 + (pB[1]-pA[1])**2)
            dist_BC = math.sqrt((pC[0]-pB[0])**2 + (pC[1]-pB[1])**2)
            dist_CA = math.sqrt((pA[0]-pC[0])**2 + (pA[1]-pC[1])**2)
            
            g_parts = []
            g_parts.append(self.extractDataAlongPath(h2, pA[0], pA[1], pB[0], pB[1], 0.0))
            g_parts.append(self.extractDataAlongPath(h2, pB[0], pB[1], pC[0], pC[1], dist_AB))
            g_parts.append(self.extractDataAlongPath(h2, pC[0], pC[1], pA[0], pA[1], dist_AB + dist_BC))
            g_parts.append(self.extractDataAlongPath(h2, pA[0], pA[1], pD[0], pD[1], dist_AB + dist_BC + dist_CA))
            
            x_buf, y_buf, ex_buf, ey_buf = [], [], [], []
            for g in g_parts:
                if not g: continue
                for k in range(g.GetN()):
                    x_buf.append(g.GetX()[k]); y_buf.append(g.GetY()[k]); ex_buf.append(0); ey_buf.append(g.GetEY()[k])
            
            if not x_buf: continue
            import array
            g_total = ROOT.TGraphErrors(len(x_buf), array.array('d', x_buf), array.array('d', y_buf), array.array('d', ex_buf), array.array('d', ey_buf))
            
            # --- Style Settings (Line) ---
            g_total.SetLineColor(ds.color)
            g_total.SetLineWidth(2)
            g_total.SetMarkerColor(ds.color)
            g_total.SetMarkerStyle(ds.markerStyle)
            
            # --- Style Settings (Error Band) ---
            # グラフを複製してエラーバンド用にする
            g_band = g_total.Clone(f"band_{ds.name}")
            g_band.SetFillColorAlpha(ds.color, 0.3) # 透過度を設定 (0.3程度が適当)
            g_band.SetFillStyle(1001) # ソリッド塗りつぶし
            g_band.SetLineWidth(0)    # バンドの枠線はなし
            
            # MultiGraphに追加 (順序: バンド -> 線)
            mg_combined.Add(g_band, "3") # "3" はエラーバンド描画
            mg_combined.Add(g_total, "L") # "L" は線のみ描画
            
            l.AddEntry(g_total, ds.name, "l") # レジェンドは線で登録
            has_data = True
            
            # 保持リストに追加
            kept_graphs.append(g_total)
            kept_bands.append(g_band)

        if not has_data: return
        mg_combined.Draw("A"); l.Draw()
        c.SaveAs(f"{self.output_dirs['Path']}/CombinedPathAnalysis_{suffix}.pdf")

if __name__ == "__main__":
    plotter = plot_Check("modeling_check_v1.json")
    try:
        plotter.run_plotCheck()
    except Exception as e:
        log.critical(f"Error: {e}", exc_info=True)
    finally:
        plotter.cleanupResources()