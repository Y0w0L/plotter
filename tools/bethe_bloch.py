import numpy as np
import matplotlib.pyplot as plt

# --- 定数 ---

# シリコンにおける電子-正孔対の生成エネルギー (eV)
ENERGY_PER_PAIR_SI = 3.6  # eV

# 各物質の物性値
# 密度 (rho) [g/cm^3]
# 質量衝突エネルギー損失 (dEdx_mass) [MeV*cm^2/g] (MIPのフェルミプラトー値)
MATERIALS = {
    'Silicon': {
        'rho': 2.33,
        'dEdx_mass': 1.664,
        'marker': '^',
        'color': 'darkorange'
    },
    'Aluminum': {
        'rho': 2.70,
        'dEdx_mass': 1.62,
        'marker': 'v',
        'color': 'red'
    },
    'Mylar': {
        'rho': 1.40,
        'dEdx_mass': 1.95,
        'marker': 'D',
        'color': 'blue'
    },
    'Scinti': { # Scintillator (Polystyrene)
        'rho': 1.05,
        'dEdx_mass': 1.94,
        'marker': 's',
        'color': 'green'
    },
    'Air': {
        'rho': 0.001225,
        'dEdx_mass': 1.83,
        'marker': 'o',
        'color': 'grey'
    }
}

# X軸の範囲: 0.01 mm から 1000 mm (10^-2 から 10^3)
thickness_mm = np.logspace(-2, 3, 100) # 100点
thickness_cm = thickness_mm * 0.1     # cmに変換

# 1行2列のサブプロットを作成
fig, ax = plt.subplots(1, 2, figsize=(16, 7))

# --- プロット 1: 厚さ vs. エネルギー損失 (ax[0] に描画) ---

for name, mat in MATERIALS.items():
    dedx_linear_mev_cm = mat['dEdx_mass'] * mat['rho']
    energy_loss_kev = dedx_linear_mev_cm * thickness_cm * 1000.0
    ax[0].plot(
        thickness_mm, 
        energy_loss_kev, 
        label=name, 
        color=mat['color'], 
        marker=mat['marker'],
        markersize=6,
        markevery=10
    )

ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].grid(True, which="both", linestyle='--', linewidth=0.5)
ax[0].set_xlabel('Thickness [mm]', fontsize=14)
ax[0].set_ylabel('Mean Energy Loss (Collision) [keV]', fontsize=14)
ax[0].set_title('MIP (3 GeV/c e-) Energy Loss vs. Thickness', fontsize=16)
ax[0].legend(title="Materials", fontsize=12)

# ★★★ 追記 (プロット1) ★★★
# 10µm = 0.01mm, 15µm = 0.015mm の位置に線を引く
ax[0].axvline(x=0.01, color='black', linestyle=':', linewidth=1.5)
ax[0].axvline(x=0.015, color='black', linestyle=':', linewidth=1.5)

# テキストラベルを追加 (y軸はlogスケールなので、下端からの倍率で位置を指定)
ymax0 = ax[0].get_ylim()[1] # y軸の最小値を取得
ax[0].text(0.01, ymax0 * 0.8, ' 10 µm', ha='left', va='top', fontsize=10,
           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))
ax[0].text(0.015, ymax0 * 0.45, ' 15 µm', ha='left', va='top', fontsize=10,
           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))
# ★★★ 追記ここまで ★★★


# --- プロット 2: 厚さ vs. 生成キャリアペア数 (ax[1] に描画) ---

mat_si = MATERIALS['Silicon']
dedx_linear_mev_cm_si = mat_si['dEdx_mass'] * mat_si['rho']
energy_loss_ev_si = dedx_linear_mev_cm_si * thickness_cm * 1e6
carrier_pairs = energy_loss_ev_si / ENERGY_PER_PAIR_SI

ax[1].plot(
    thickness_mm, 
    carrier_pairs, 
    label='Silicon',
    color=mat_si['color'],
    marker=mat_si['marker'],
    markersize=6,
    markevery=10
)

ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].grid(True, which="both", linestyle='--', linewidth=0.5)
ax[1].set_xlabel('Thickness [mm]', fontsize=14)
ax[1].set_ylabel('Mean Carrier Pairs (e-h pairs)', fontsize=14)
ax[1].set_title('MIP (3 GeV/c e-) Carrier Generation in Silicon', fontsize=16)
ax[1].legend(fontsize=12)

# ★★★ 追記 (プロット2) ★★★
ax[1].axvline(x=0.01, color='black', linestyle=':', linewidth=1.5)
ax[1].axvline(x=0.015, color='black', linestyle=':', linewidth=1.5)

# テキストラベルを追加
ymax1 = ax[1].get_ylim()[1] # y軸の最小値を取得
ax[1].text(0.01, ymax1 * 0.8, ' 10 µm', ha='left', va='top', fontsize=10,
           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))
ax[1].text(0.015, ymax1 * 0.5, ' 15 µm', ha='left', va='top', fontsize=10,
           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))
# ★★★ 追記ここまで ★★★


# --- 最後にレイアウトを自動調整して表示 ---
plt.tight_layout()
plt.show()
plt.savefig("bethe_bloch.pdf")