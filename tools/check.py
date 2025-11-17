import ROOT
import argparse
import os
import sys

def setup_style():
    """ROOTの描画スタイルを設定する"""
    ROOT.gROOT.SetBatch(True)  # ウィンドウを表示しないバッチモード
    ROOT.gStyle.SetOptStat(1111)  # 基本統計量の表示
    ROOT.gStyle.SetOptFit(1111)   # フィッティング結果の表示 (確率, カイ二乗, エラー, パラメータ)
    ROOT.gStyle.SetStatX(0.9)     # 統計ボックスの位置調整
    ROOT.gStyle.SetStatY(0.9)
    ROOT.gStyle.SetStatW(0.2)     # 統計ボックスの幅
    ROOT.gStyle.SetStatH(0.2)     # 統計ボックスの高さ
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "")

def fit_and_draw(canvas, hist, pdf_filename, title_suffix=""):
    """ヒストグラムをガウス分布でフィットしてPDFにページを追加する"""
    if not hist or hist.GetEntries() < 10:
        # エントリが少なすぎる場合はスキップ
        return

    # キャンバスをクリア
    canvas.Clear()
    
    # タイトルなどを調整
    hist.SetTitle(f"{hist.GetTitle()}{title_suffix}")
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleOffset(1.4)
    
    # ガウス分布でフィッティング ("Q"は静音モード, "S"は結果取得)
    # 範囲をヒストグラムの表示範囲全体にするか、少し狭めるか調整可能
    # ここでは単純な "gaus" を使用
    fit_result = hist.Fit("gaus", "QS")
    
    # 描画
    hist.Draw()
    
    # フィット結果が有効なら、フィット線も描画されている
    
    # PDFに現在のキャンバスの状態を出力
    canvas.Print(pdf_filename)

def main():
    parser = argparse.ArgumentParser(description="Fit residuals with Gaussian and save to PDF.")
    parser.add_argument("-i", "--input", type=str, default="analysis_py.root", help="Input ROOT file from analysis.")
    parser.add_argument("-o", "--output", type=str, default="residuals_fitted.pdf", help="Output PDF file name.")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    setup_style()

    # ROOTファイルを開く
    tfile = ROOT.TFile.Open(args.input, "READ")
    if not tfile or tfile.IsZombie():
        print(f"Error: Cannot open ROOT file '{args.input}'.")
        sys.exit(1)

    # キャンバスの作成
    cv = ROOT.TCanvas("cv", "Residual Fits", 800, 600)
    
    # PDFの開始（ファイル名の末尾に '[' を付けると、ページを開くだけで保存しない）
    cv.Print(args.output + "[")

    print(f"Processing file: {args.input}")
    print(f"Outputting to: {args.output}")

    # クラスターサイズ 1 から 10 までのループ
    max_size = 10
    
    # 処理する軸 (Residual Rはガウス分布ではないため、通常XとYのみフィットします)
    axes = ["x", "y"] 

    for i in range(1, max_size + 2):
        # サイズ文字列の決定 (1, 2, ..., 10, 11_plus)
        size_str = str(i) if i <= max_size else f"{max_size + 1}_plus"
        
        for axis in axes:
            hist_name = f"residual_{axis}_size_{size_str}"
            hist = tfile.Get(hist_name)

            if hist:
                print(f"Fitting {hist_name}...")
                fit_and_draw(cv, hist, args.output)
            else:
                # ヒストグラムが見つからない場合（特定のサイズが存在しないなど）
                # print(f"Warning: Histogram '{hist_name}' not found.")
                pass

    # Residual R もプロットだけしたい（フィットなし）場合は以下のコメントを解除
    """
    for i in range(1, max_size + 2):
        size_str = str(i) if i <= max_size else f"{max_size + 1}_plus"
        hist_name = f"residual_r_size_{size_str}"
        hist = tfile.Get(hist_name)
        if hist and hist.GetEntries() > 0:
             cv.Clear()
             hist.SetTitle(f"{hist.GetTitle()} (No Fit)")
             hist.Draw()
             cv.Print(args.output)
    """

    # PDFの終了（ファイル名の末尾に ']' を付けると、ファイルを閉じる）
    cv.Print(args.output + "]")
    
    tfile.Close()
    print("Done.")

if __name__ == "__main__":
    main()