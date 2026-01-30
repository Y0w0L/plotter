import ROOT
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def get_path_points(pitch):
    """Calculate cumulative distance of the path based on pixel pitch"""
    hp = pitch / 2.0
    l_ab, l_bc, l_ca, l_ad = hp, hp, np.sqrt(hp**2 + hp**2), hp
    points = [0, l_ab, l_ab + l_bc, l_ab + l_bc + l_ca, l_ab + l_bc + l_ca + l_ad]
    labels = ["A", "B", "C", "A", "D"]
    return points, labels

def draw_path_inset(ax):
    """Draw a schematic diagram of the path in the upper right corner"""
    ax_ins = inset_axes(ax, width="22%", height="22%", loc='upper right', 
                        bbox_to_anchor=(0, -0.18, 1, 1), bbox_transform=ax.transAxes)
    
    rect = patches.Rectangle((-1.0, -1.0), 2.0, 2.0, linewidth=1, 
                             edgecolor='gray', facecolor='none', linestyle=':', alpha=0.7)
    ax_ins.add_patch(rect)

    A, B, C, D = (0, 0), (0, 1.0), (1.0, 1.0), (1.0, 0)
    arrow_props = dict(arrowstyle='->', lw=2, mutation_scale=12)
    path_colors = ["#ff7f0e", "#17becf", "#d62728", "#8c564b"] 
    
    ax_ins.annotate('', xy=B, xytext=A, arrowprops=dict(**arrow_props, color=path_colors[0]))
    ax_ins.annotate('', xy=C, xytext=B, arrowprops=dict(**arrow_props, color=path_colors[1]))
    ax_ins.annotate('', xy=A, xytext=C, arrowprops=dict(**arrow_props, color=path_colors[2]))
    ax_ins.annotate('', xy=D, xytext=A, arrowprops=dict(**arrow_props, color=path_colors[3]))

    label_font = {'fontsize': 10, 'fontweight': 'bold'}
    ax_ins.text(0, -0.15, 'A', **label_font, ha='center', va='top')
    ax_ins.text(0, 1.1, 'B', **label_font, ha='center', va='bottom')
    ax_ins.text(1.1, 1.1, 'C', **label_font, ha='left', va='bottom')
    ax_ins.text(1.1, 0, ' D', **label_font, ha='left', va='center')

    ax_ins.set_xlim(-1.3, 1.5); ax_ins.set_ylim(-1.3, 1.5)
    ax_ins.set_aspect('equal'); ax_ins.axis('off')

def compare_in_pixel_paths(configs):
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Red system colors (SPS) and Blue system colors (SIM)
    sps_colors = ["#fc1515", "#ff0000", "#cc0000", "#990000"]
    sim_colors = ["#0981f8", "#0066ff", "#0047b3", "#002966"]

    target_pitch = configs[0]["pitch"]
    points, labels = get_path_points(target_pitch)

    # 1. Vertical guide lines and secondary X-axis for labels
    for p in points:
        ax.axvline(x=p, color="gray", linestyle="--", alpha=0.3)
    
    secax = ax.secondary_xaxis('top')
    secax.set_xticks(points)
    secax.set_xticklabels(labels, fontweight='bold', fontsize=12, color='gray')
    secax.tick_params(axis='x', direction='in', length=0)

    # 2. Data Plotting
    for conf in configs:
        f = ROOT.TFile.Open(conf["path"])
        if not f or f.IsZombie(): 
            print(f"Error: Could not open {conf['path']}")
            continue
        
        mg = f.Get(f"{conf['qty']}/{conf['cond']}/mg_combined_path_{conf['qty']}")
        if not mg: 
            print(f"Error: MultiGraph not found in {conf['path']}")
            f.Close(); continue
            
        graphs = mg.GetListOfGraphs()
        colors = sps_colors if conf["marker"] == "^" else sim_colors

        x_all, y_all, yerr_all = [], [], []
        for g in graphs:
            n = g.GetN()
            if n == 0: continue
            x_all.extend(np.ndarray(n, 'd', g.GetX()))
            y_all.extend(np.ndarray(n, 'd', g.GetY()))
            yerr_all.extend(np.ndarray(n, 'd', g.GetEY()) if hasattr(g, 'GetEY') else np.zeros(n))
        
        x_all, y_all, yerr_all = np.array(x_all), np.array(y_all), np.array(yerr_all)

        for i in range(len(points) - 1):
            mask = (x_all >= points[i] - 1e-6) & (x_all <= points[i+1] + 1e-6)
            if not np.any(mask): continue
            sort_idx = np.argsort(x_all[mask])
            ax.errorbar(x_all[mask][sort_idx], y_all[mask][sort_idx], 
                        yerr=yerr_all[mask][sort_idx], fmt=conf["marker"], 
                        color=colors[i % 4], markersize=9, capsize=3, 
                        linestyle='-', alpha=0.8, label="_nolegend_")
            ax.fill_between(x_all[mask][sort_idx], 
                            y_all[mask][sort_idx] - yerr_all[mask][sort_idx], 
                            y_all[mask][sort_idx] + yerr_all[mask][sort_idx], 
                            color=colors[i % 4], alpha=0.3, linewidth=0, label="_nolegend_")

        f.Close()

    # 3. Legend, Chip Info and Inset
    # Construct chip information text from configs
    info_text = ""
    for conf in configs:
        prefix = "SPS" if conf["marker"] == "^" else "SIM"
        info_text += f"{prefix} ({conf.get('chip', 'N/A')}): {conf['cond']}\n"
    
    # Place text box in the upper left
    ax.text(0.02, 0.95, info_text.strip(), transform=ax.transAxes, 
            fontsize=12, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))

    chip_raw = configs[0].get("chip", "N/A")
    chip_label = chip_raw.replace("p22p5", "p22.5").replace("_", " ").replace("gap", "GAP").replace("std", "STD")

    legend_elements = [
        Line2D([0], [0], color='gray', marker='^', linestyle='-', label=f'SPS {chip_label}'),
        Line2D([0], [0], color='gray', marker='o', linestyle='-', label=f'SIM {chip_label}')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, frameon=True)
    draw_path_inset(ax)

    # 4. Axis Styling
    ax.tick_params(direction='in', top=True, right=True, which='both', length=6, labelsize=16)
    ax.set_xlabel("Distance along path [$\mu$m]", fontsize=18, labelpad=12)
    if configs[0]["qty"] == "residual":
        y_label = r"$r_{\mathrm{exp}} - r_{\mathrm{hit}} \ [\mu\mathrm{m}]$"
    else:
        y_label = configs[0]["qty"].replace("_", " ").title()
    ax.set_ylabel(y_label, fontsize=18, labelpad=12)
    ax.grid(True, linestyle=':', alpha=0.5)
    
    plt.subplots_adjust(top=0.95, bottom=0.13, left=0.1, right=0.97)
    
    chip_name = configs[0]["chip"]
    output_pdf = f"path_{qty}_{chip_name}.pdf"
    plt.savefig(output_pdf)
    print(f"Saved: {output_pdf}")
    plt.close() # Close to prevent memory issues during batch processing

if __name__ == "__main__":
    # プロットしたい項目のリスト
    quantities = ["cluster_size", "residual"]
    # 比較したいチップ構成のリスト
    chips = ["p15_gap", "p15_std", "p22p5_gap", "p22p5_std"]

    for qty in quantities:
        for chip in chips:
            # チップ名からピッチを判定 (p15 -> 15.0, p22p5 -> 22.5)
            pitch = 15.0 if "p15" in chip else 22.5
            
            # 各チップに対応するファイルパスを動的に生成
            configs = [
                {
                    "path": f"/home/towa/alice3/plotter/plot/inPixelPathAnalysis_sps202404_{chip}.root",
                    "qty": qty,
                    "cond": "Seed1000_Neighbor600",
                    "chip": chip,
                    "pitch": pitch,
                    "marker": "^"
                },
                {
                    "path": f"/home/towa/alice3/plotter/plot/inPixelPathAnalysis_ce65inPixel2601_{chip}.root",
                    "qty": qty,
                    "cond": "Seed240_Neighbor144",
                    "chip": chip,
                    "pitch": pitch,
                    "marker": "o"
                }
            ]
            
            print(f"Processing: {qty} for {chip}...")
            compare_in_pixel_paths(configs)