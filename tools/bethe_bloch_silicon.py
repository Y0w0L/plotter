import sys
import numpy as np # (計算範囲を広げるために numpy をインポートしておくと便利です)
import pandas as pd # (データフレーム化やCSV保存用にインポート)
import matplotlib.pyplot as plt # (プロット用にインポート)
import matplotlib.ticker as ticker # (プロットの目盛り調整用にインポート)


# --- 定数 (シリコン用) ---

# シリコンの密度 (g/cm^3)
RHO_SI = 2.33  # g/cm^3

# 最小電離粒子 (MIP) の平均衝突エネルギー損失 (dE/dx)
# (ベーテ・ブロッホの式によるフェルミプラトーの値)
# 120GeV/c の π+ も 3GeV/c の e- も、MIPとして扱えます。
DEDX_COLL_MIP = 1.664  # MeV * cm^2 / g

# 1ペアの電子・ホール対を生成するのに必要なエネルギー (eV)
ENERGY_PER_PAIR_SI = 3.6  # eV

def calculate_energy_loss(thickness_microns: float) -> dict:
    """
    シリコン検出器中でのMIP（120GeV/c π+）の平均衝突エネルギー損失を計算する。
    
    Args:
        thickness_microns (float): シリコンの厚さ (マイクロメートル, µm)
        
    Returns:
        dict: 計算結果を含む辞書
    """
    if thickness_microns <= 0:
        return {"error": "厚さは正の値である必要があります"}

    # --- 単位換算 ---
    # µm -> cm
    thickness_cm = thickness_microns * 1e-4  # 1 µm = 1e-4 cm
    
    # --- エネルギー損失 (ΔE) の計算 ---
    # dE/dx (MeV/cm) を計算
    # dE/dx [MeV/cm] = dE/dx [MeV*cm^2/g] * ρ [g/cm^3]
    dedx_per_cm = DEDX_COLL_MIP * RHO_SI  # 約 3.88 MeV/cm
    
    # 平均エネルギー損失 ΔE (MeV)
    # ΔE [MeV] = dE/dx [MeV/cm] * Δx [cm]
    delta_e_mev = dedx_per_cm * thickness_cm
    
    # MeV を keV に変換
    delta_e_kev = delta_e_mev * 1000.0
    
    # --- 生成電荷の計算 ---
    # ΔE (eV) に変換
    delta_e_ev = delta_e_mev * 1e6
    
    # 生成される電子・ホール対の平均数
    avg_pairs = delta_e_ev / ENERGY_PER_PAIR_SI
    
    return {
        "thickness_microns": thickness_microns,
        "thickness_cm": thickness_cm,
        "dedx_mev_per_cm": dedx_per_cm,
        "mean_energy_loss_mev": delta_e_mev,
        "mean_energy_loss_kev": delta_e_kev,
        "avg_electron_hole_pairs": avg_pairs
    }

# --- メインの実行部分 ---
if __name__ == "__main__":
    
    # デフォルトで計算したい厚さ (µm)
    # (前回のリクエストに基づき 5um から 20um まで 0.5um ステップで計算する例)
    thicknesses_to_test = np.arange(5.0, 50.5, 0.5) 
    
    # (オプション) コマンドラインから引数として厚さを渡すこともできます
    # 例: python calculate_loss.py 50 100
    if len(sys.argv) > 1:
        try:
            thicknesses_to_test = [float(arg) for arg in sys.argv[1:]]
            print(f"コマンドライン引数で指定された厚さで計算します: {thicknesses_to_test}")
        except ValueError:
            print(f"使用法: python {sys.argv[0]} [厚さ1] [厚さ2] ...")
            print(f"例: python {sys.argv[0]} 50 300")
            sys.exit(1)

    print("--- 120GeV/c π+ (MIP) のシリコン中での衝突エネルギー損失計算 ---")
    
    results_list = []
    
    for thickness in thicknesses_to_test:
        result = calculate_energy_loss(thickness)
        
        if "error" in result:
            print(f"\n[厚さ: {thickness} µm] エラー: {result['error']}")
        else:
            results_list.append(result)
            # コンソールへの出力はオプション (データが多い場合はコメントアウト推奨)
            # print(f"\n--- [厚さ: {thickness:.1f} µm] ---")
            # print(f"  平均エネルギー損失 (ΔE): {result['mean_energy_loss_kev']:.2f} keV")
            # print(f"  平均生成電荷 (e-hペア): {result['avg_electron_hole_pairs']:,.0f} ペア")

    # --- DataFrame化とプロット (前回のコードから流用) ---
    if not results_list:
        print("計算結果がありません。")
        sys.exit(1)
        
    df = pd.DataFrame(results_list)
    
    print("\n--- 計算結果 (先頭5件) ---")
    print(df.head())
    print("...")
    print("\n--- 計算結果 (末尾5件) ---")
    print(df.tail())
    
    # 結果をCSVファイルに保存
    # csv_filename = "silicon_energy_loss_120GeV_pion_MIP.csv"
    # df.to_csv(csv_filename, index=False, float_format='%.4f')
    # print(f"\n計算結果全体を {csv_filename} に保存しました。")

    # --- プロットの作成 ---
    plt.figure(figsize=(10, 6))
    plt.plot(df['thickness_microns'], df['mean_energy_loss_kev'], marker='.')
    plt.xlabel("Silicon Thickness (µm)")
    plt.ylabel("Mean Energy Loss (keV)")
    plt.title("Mean Energy Loss (120 GeV/c $\pi^+$ MIP) vs. Silicon Thickness")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    # x軸の目盛りを調整 (入力範囲に応じて調整が必要な場合があります)
    if max(thicknesses_to_test) - min(thicknesses_to_test) <= 20:
         plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(2))
         plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    
    plt.tight_layout()
    plot_filename_kev = "energy_loss_vs_thickness_120GeV_pion.pdf"
    plt.savefig(plot_filename_kev)
    print(f"エネルギー損失のグラフを {plot_filename_kev} に保存しました。")
    plt.clf() # グラフをクリア

    # 2. 厚さ vs 平均生成電荷 (e-hペア)
    plt.figure(figsize=(10, 6))
    plt.plot(df['thickness_microns'], df['avg_electron_hole_pairs'], marker='.', color='orange')
    plt.xlabel("Silicon Thickness (µm)")
    plt.ylabel("Average e-h Pairs")
    plt.title("Average e-h Pairs (120 GeV/c $\pi^+$ MIP) vs. Silicon Thickness")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    if max(thicknesses_to_test) - min(thicknesses_to_test) <= 20:
        plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(2))
        plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    plt.tight_layout()
    plot_filename_pairs = "pairs_vs_thickness_120GeV_pion.pdf"
    plt.savefig(plot_filename_pairs)
    print(f"生成電荷数のグラフを {plot_filename_pairs} に保存しました。")
    plt.clf() # グラフをクリア