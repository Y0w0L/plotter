import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# ==========================================
# 1. Recombination / Trapping Models
# ==========================================

def srh_model(N, temp_k, particle_type):
    if particle_type == 'electron':
        tau_0 = 1.0e-5
        Nd0 = 1.0e16
    elif particle_type == 'hole':
        tau_0 = 4.0e-4
        Nd0 = 7.1e15
    else:
        return np.zeros_like(N)
    tau_N = tau_0 / (1 + N / Nd0)
    tau_final = tau_N * (300.0 / temp_k)**1.5
    return tau_final

def auger_model(N, particle_type):
    if particle_type == 'electron':
        Ca = 2.8e-31
    elif particle_type == 'hole':
        Ca = 9.9e-32
    else:
        return np.zeros_like(N)
    tau = 1.0 / (Ca * N**2)
    return tau

def combined_recombination_model(N, temp_k, particle_type):
    tau_srh = srh_model(N, temp_k, particle_type)
    tau_auger = auger_model(N, particle_type)
    inv_tau = (1.0 / tau_srh) + (1.0 / tau_auger)
    return 1.0 / inv_tau

def trapping_models(fluence, temp_k, model_name, particle_type):
    if model_name == 'Ljubljana':
        T0 = 263.15
        if particle_type == 'electron': params = (5.6e-16, -0.86)
        else: params = (7.7e-16, -1.52)
        beta = params[0] * (temp_k/T0)**params[1]
        return 1.0 / (beta * fluence)
    elif model_name == 'Dortmund':
        if particle_type == 'electron': gamma = 5.13e-16
        else: gamma = 5.04e-16
        return 1.0 / (gamma * fluence)
    elif model_name == 'CMS Tracker':
        if particle_type == 'electron': p = (1.71e-16, -0.114)
        else: p = (2.79e-16, -0.093)
        inv_tau = p[0] * fluence + p[1]
        return 1.0 / inv_tau
    elif model_name == 'Mandic':
        if particle_type == 'electron': c, k = (0.54, -0.62)
        else: c, k = (0.0427, -0.62)
        return c * np.power(fluence, k)

# ==========================================
# 2. Mobility Models
# ==========================================

def jacoboni_canali_model(E, temp_k, particle_type):
    # E: Electric Field [V/cm]
    if particle_type == 'electron':
        vm = 1.53e9 * temp_k**(-0.87)
        Ec = 1.01 * temp_k**1.55
        beta = 2.57e-2 * temp_k**0.66
    else: # hole
        vm = 1.62e8 * temp_k**(-0.52)
        Ec = 1.24 * temp_k**1.68
        beta = 0.46 * temp_k**0.17
    
    mu = (vm / Ec) / (1 + (E / Ec)**beta)**(1/beta)
    return mu

def canali_model(E, temp_k, particle_type):
    # Parameters similar to Jacoboni but different vm for electrons
    if particle_type == 'electron':
        vm = 1.43e9 * temp_k**(-0.87) # Canali specific
        Ec = 1.01 * temp_k**1.55      # Same as Jacoboni
        beta = 2.57e-2 * temp_k**0.66 # Same as Jacoboni
    else: # hole (Same as Jacoboni)
        vm = 1.62e8 * temp_k**(-0.52)
        Ec = 1.24 * temp_k**1.68
        beta = 0.46 * temp_k**0.17

    mu = (vm / Ec) / (1 + (E / Ec)**beta)**(1/beta)
    return mu, vm, beta # Return vm, beta for Extended Canali usage

def masetti_model(N, temp_k, particle_type):
    # Calculates low-field mobility based on doping
    if particle_type == 'electron':
        mu0 = 68.5
        mumax = 1414 * (temp_k / 300.0)**(-2.5)
        Cr = 9.20e16
        alpha = 0.711
        mu1 = 56.1
        Cs = 3.41e20
        beta_p = 1.98
        
        term1 = (mumax - mu0) / (1 + (N / Cr)**alpha)
        term2 = mu1 / (1 + (Cs / N)**beta_p)
        return mu0 + term1 - term2
    else: # hole
        mu0 = 44.9
        mumax = 470.5 * (temp_k / 300.0)**(-2.2)
        Cr = 2.23e17
        alpha = 0.719
        mu1 = 29.0
        Cs = 6.1e20
        beta_p = 2.0
        Pc = 9.23e16
        
        term1 = mu0 * np.exp(-Pc / N)
        term2 = mumax / (1 + (N / Cr)**alpha)
        term3 = mu1 / (1 + (Cs / N)**beta_p)
        return term1 + term2 - term3

def arora_model(N, temp_k, particle_type):
    # Calculates low-field mobility based on doping
    if particle_type == 'electron':
        mumin = 88.0 * (temp_k / 300.0)**(-0.57)
        mu0_p = 7.40e8 * temp_k**(-2.33)
        Nref = 1.26e17 * (temp_k / 300.0)**2.4
    else: # hole
        mumin = 54.3 * (temp_k / 300.0)**(-0.57)
        mu0_p = 1.36e8 * temp_k**(-2.23)
        Nref = 2.35e17 * (temp_k / 300.0)**2.4
    
    alpha = 0.88 * (temp_k / 300.0)**(-0.146)
    return mumin + mu0_p / (1 + (N / Nref)**alpha)

def extended_canali_model(E, N, temp_k, particle_type):
    # Uses Masetti for low field, Canali parameters for saturation
    mu_low = masetti_model(N, temp_k, particle_type)
    _, vm, beta = canali_model(1.0, temp_k, particle_type) # get params
    
    mu = mu_low / (1 + (mu_low * E / vm)**beta)**(1/beta)
    return mu

def arora_canali_model(E, N, temp_k, particle_type):
    # Uses Arora for low field, Canali parameters for saturation
    mu_low = arora_model(N, temp_k, particle_type)
    _, vm, beta = canali_model(1.0, temp_k, particle_type) # get params
    
    mu = mu_low / (1 + (mu_low * E / vm)**beta)**(1/beta)
    return mu

def hamburg_model(E, temp_k, particle_type, high_field_version=False):
    # Hamburg and Hamburg High-Field models share similar structure
    if particle_type == 'electron':
        if not high_field_version:
            mu0 = 1530 * (temp_k / 300.0)**(-2.42)
            vsat = 1.03e7 * (temp_k / 300.0)**(-0.226)
        else: # High Field
            mu0 = 1430 * (temp_k / 300.0)**(-1.99)
            vsat = 1.05e7 * (temp_k / 300.0)**(-0.302)
        
        # 1/mu = 1/mu0 + E/vsat -> mu = 1 / (1/mu0 + E/vsat)
        inv_mu = 1.0/mu0 + E/vsat
        return 1.0 / inv_mu

    else: # hole
        if not high_field_version:
            mu0 = 464 * (temp_k / 300.0)**(-2.20)
            b = 9.57e-8 * (temp_k / 300.0)**(-0.101)
            c = -3.31e-13
            E0 = 2640 * (temp_k / 300.0)**0.526
        else: # High Field
            mu0 = 457 * (temp_k / 300.0)**(-2.80)
            b = 9.57e-8 * (temp_k / 300.0)**(-0.155)
            c = -3.24e-13
            E0 = 2970 * (temp_k / 300.0)**0.563

        # Vectorized calculation for E array
        mu = np.zeros_like(E)
        mask_low = E < E0
        mask_high = E >= E0
        
        mu[mask_low] = mu0
        
        # For E >= E0: 1/mu = 1/mu0 + b(E-E0) + c(E-E0)^2
        term_high = 1.0/mu0 + b * (E[mask_high] - E0) + c * (E[mask_high] - E0)**2
        mu[mask_high] = 1.0 / term_high
        
        return mu

# ==========================================
# 3. Simulation Setup
# ==========================================

temp_c = 30.0
temp_k = 273.15 + temp_c

# --- Data for Recombination Plot ---
N_doping = np.logspace(14, 20, 100)
srh_e = srh_model(N_doping, temp_k, 'electron')
srh_h = srh_model(N_doping, temp_k, 'hole')
aug_e = auger_model(N_doping, 'electron')
aug_h = auger_model(N_doping, 'hole')
comb_e = combined_recombination_model(N_doping, temp_k, 'electron')
comb_h = combined_recombination_model(N_doping, temp_k, 'hole')

# --- Data for Trapping Plot ---
fluence = np.linspace(1e12, 1e16, 1000)
trap_data = []
trap_models_list = ['Ljubljana', 'Dortmund', 'CMS Tracker', 'Mandic']
for m in trap_models_list:
    trap_data.append((m, 'Electron', trapping_models(fluence, temp_k, m, 'electron')))
    trap_data.append((m, 'Hole', trapping_models(fluence, temp_k, m, 'hole')))

# --- Data for Mobility/Velocity Plot ---
E_field = np.linspace(0, 10000, 500) # 0 to 10 kV/cm
N_fixed = 1e14 # Fixed doping concentration for Extended/Arora models [cm^-3]

# Calculate Mobilities
mob_jac_e = jacoboni_canali_model(E_field, temp_k, 'electron')
mob_jac_h = jacoboni_canali_model(E_field, temp_k, 'hole')

mob_can_e, _, _ = canali_model(E_field, temp_k, 'electron') # Canali returns tuple
# Note: Canali Hole is identical to Jacoboni Hole in this formulation, so typically not plotted separately or overlaps

mob_ext_e = extended_canali_model(E_field, N_fixed, temp_k, 'electron')
mob_ext_h = extended_canali_model(E_field, N_fixed, temp_k, 'hole')

mob_arora_e = arora_canali_model(E_field, N_fixed, temp_k, 'electron')
mob_arora_h = arora_canali_model(E_field, N_fixed, temp_k, 'hole')

mob_ham_e = hamburg_model(E_field, temp_k, 'electron', high_field_version=False)
mob_ham_h = hamburg_model(E_field, temp_k, 'hole', high_field_version=False)

mob_ham_hf_e = hamburg_model(E_field, temp_k, 'electron', high_field_version=True)
mob_ham_hf_h = hamburg_model(E_field, temp_k, 'hole', high_field_version=True)

# Common Plot Settings
plt.rcParams['font.family'] = 'sans-serif'

# ==========================================
# 4. Plot Generation
# ==========================================

# --- Plot 1: Recombination Lifetime ---
fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(N_doping, comb_e, label='Combined (Electron)', color='navy', linestyle='-', linewidth=2)
ax1.plot(N_doping, srh_e, label='SRH (Electron)', color='tab:blue', linestyle='--')
ax1.plot(N_doping, aug_e, label='Auger (Electron)', color='dodgerblue', linestyle=':')
ax1.plot(N_doping, comb_h, label='Combined (Hole)', color='darkred', linestyle='-', linewidth=2)
ax1.plot(N_doping, srh_h, label='SRH (Hole)', color='tab:red', linestyle='--')
ax1.plot(N_doping, aug_h, label='Auger (Hole)', color='tomato', linestyle=':')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'Doping Concentration $N$ [cm$^{-3}$]', fontsize=12)
ax1.set_ylabel(r'Carrier Lifetime $\tau$ [s]', fontsize=12)
ax1.grid(True, linestyle='--', alpha=0.7, which='both')
ax1.set_xlim(1e14, 1e20)
ax1.set_ylim(1e-10, 1e-2) 
ax1.text(0.02, 0.95, "Recombination Lifetime vs. Doping", transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
ax1.legend(loc='lower left', framealpha=0.9, fontsize=9)
plt.tight_layout()
fig1.savefig("recombination_lifetime.pdf")
plt.close(fig1)

# --- Plot 2: Trapping Time ---
fig2, ax2 = plt.subplots(figsize=(8, 6))
model_colors = {'Ljubljana': 'tab:blue', 'Dortmund': 'tab:purple', 'CMS Tracker': 'tab:red', 'Mandic': 'tab:green'}
for name, p_type, data in trap_data:
    c = model_colors[name]
    s = '-' if p_type == 'Electron' else '--'
    ax2.plot(fluence, data, label=f"{name} ({p_type})", color=c, linestyle=s)
ax2.set_yscale('log')
ax2.set_xlabel(r'Fluence $\Phi_{eq}$ [n$_{eq}$/cm$^2$]', fontsize=12)
ax2.set_ylabel(r'Trapping Time $\tau$ [ns]', fontsize=12)
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax2.set_ylim(1e-12, 1e5)
ax2.text(0.02, 0.95, "Trapping Time vs. Fluence", transform=ax2.transAxes, fontsize=14, fontweight='bold', va='center')
ax2.legend(loc='upper right', title='Radiation Damage Models', framealpha=0.9, fontsize=9, ncol=2)
plt.tight_layout()
fig2.savefig("trapping_time.pdf")
plt.close(fig2)

# --- Plot 3: Mobility and Velocity vs Electric Field ---
fig3, (ax_mob, ax_vel) = plt.subplots(1, 2, figsize=(16, 6))

# --- Left: Mobility vs E-field ---
# Electron Models (Cold Colors: Blue, Cyan, Green, Purple)
ax_mob.plot(E_field, mob_jac_e, label='Electron/Jacoboni-Canali', color='blue', linestyle='-')
ax_mob.plot(E_field, mob_can_e, label='Electron/Canali', color='cyan', linestyle='-')
ax_mob.plot(E_field, mob_ext_e, label='Electron/Extended Canali', color='forestgreen', linestyle='-')
ax_mob.plot(E_field, mob_arora_e, label='Electron/Arora-Canali', color='limegreen', linestyle='-')
ax_mob.plot(E_field, mob_ham_e, label='Electron/Hamburg', color='darkviolet', linestyle='-')
ax_mob.plot(E_field, mob_ham_hf_e, label='Electron/Hamburg High-Field', color='magenta', linestyle='-')

# Hole Models (Warm Colors: Red, Orange, Gold, Pink)
ax_mob.plot(E_field, mob_jac_h, label='Hole/Jacoboni-Canali', color='red', linestyle='-')
# Canali hole is identical to Jacoboni hole
ax_mob.plot(E_field, mob_ext_h, label='Hole/Extended Canali', color='darkorange', linestyle='-')
ax_mob.plot(E_field, mob_arora_h, label='Hole/Arora-Canali', color='gold', linestyle='-')
ax_mob.plot(E_field, mob_ham_h, label='Hole/Hamburg', color='deeppink', linestyle='-')
ax_mob.plot(E_field, mob_ham_hf_h, label='Hole/Hamburg High-Field', color='hotpink', linestyle='-')

ax_mob.set_xlabel('Electric Field [V/cm]', fontsize=12)
ax_mob.set_ylabel(r'Charge Carrier Mobility [cm$^2$/Vs]', fontsize=12)
ax_mob.grid(True, linestyle='--', alpha=0.7)
ax_mob.text(0.02, 0.95, "Charge Carrier Mobility vs. Electric Field", transform=ax_mob.transAxes, fontsize=14, fontweight='bold')
ax_mob.text(0.02, 0.90, f"Plotted on {datetime.now().strftime('%d %b %Y')}", transform=ax_mob.transAxes, fontsize=10)
ax_mob.text(0.02, 0.86, f"N = {N_fixed:.0e} [cm$^{{-3}}$]", transform=ax_mob.transAxes, fontsize=10)
ax_mob.legend(title='Models', fontsize=8, loc='center right')


# --- Right: Velocity vs E-field (v = mu * E) ---
# Electron Models
ax_vel.plot(E_field, mob_jac_e * E_field, label='Electron/Jacoboni-Canali', color='tab:blue')
ax_vel.plot(E_field, mob_can_e * E_field, label='Electron/Canali', color='tab:green')
ax_vel.plot(E_field, mob_ext_e * E_field, label='Electron/Extended Canali', color='tab:red')
ax_vel.plot(E_field, mob_arora_e * E_field, label='Electron/Arora-Canali', color='tab:cyan')
ax_vel.plot(E_field, mob_ham_e * E_field, label='Electron/Hamburg', color='tab:purple')
ax_vel.plot(E_field, mob_ham_hf_e * E_field, label='Electron/Hamburg High-Field', color='tab:brown')

# Hole Models
ax_vel.plot(E_field, mob_jac_h * E_field, label='Hole/Jacoboni-Canali', color='tab:orange')
ax_vel.plot(E_field, mob_ext_h * E_field, label='Hole/Extended Canali', color='tab:pink')
ax_vel.plot(E_field, mob_arora_h * E_field, label='Hole/Arora-Canali', color='tab:olive')
ax_vel.plot(E_field, mob_ham_h * E_field, label='Hole/Hamburg', color='tab:gray')
ax_vel.plot(E_field, mob_ham_hf_h * E_field, label='Hole/Hamburg High-Field', color='gold')

ax_vel.set_xlabel('Electric Field [V/cm]', fontsize=12)
ax_vel.set_ylabel('Carrier Velocity [cm/s]', fontsize=12)
ax_vel.grid(True, linestyle='--', alpha=0.7)
ax_vel.text(0.02, 0.95, "Charge Carrier Velocity vs. Electric Field", transform=ax_vel.transAxes, fontsize=14, fontweight='bold')
ax_vel.text(0.02, 0.90, f"Plotted on {datetime.now().strftime('%d %b %Y')}", transform=ax_vel.transAxes, fontsize=10)
ax_vel.legend(title='Models', fontsize=8, loc='lower right')
ax_vel.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.tight_layout()
fig3.savefig("mobility_velocity_models.pdf")
plt.close(fig3)

print("Saved: recombination_lifetime.pdf, trapping_time.pdf, and mobility_velocity_models.pdf")