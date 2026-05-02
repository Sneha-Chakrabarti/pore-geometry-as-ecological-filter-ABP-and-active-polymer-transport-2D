import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from scipy.stats import expon

# make sure src/ is on path
sys.path.insert(0, os.path.dirname(__file__))
from abp_porous import (PorousMedia, ABPSimulation,
                        compute_msd, classify_frames,
                        trapping_probability, effective_diffusivity,
                        crossover_length)
from active_polymer import (ActivePolymer, run_msd_sweep, run_rg_sweep)

# ─────────────────────────────────────────────
# Global aesthetic settings
# ─────────────────────────────────────────────
PALETTE = {
    "silt"  : "#E07B54",
    "loam"  : "#5B8DB8",
    "sand"  : "#4CAF78",
    "linear": "#D4608A",
    "ring"  : "#7B62C7",
    "flex"  : "#3AAFB9",
    "semi"  : "#F4A261",
}

plt.rcParams.update({
    "figure.dpi"       : 150,
    "font.family"      : "serif",
    "font.size"        : 11,
    "axes.spines.top"  : False,
    "axes.spines.right": False,
    "axes.labelsize"   : 12,
    "legend.frameon"   : False,
    "lines.linewidth"  : 2.0,
})

RESULTS = "results"
os.makedirs(RESULTS, exist_ok=True)

FAST = "--fast" in sys.argv

# ─────────────────────────────────────────────
# Simulation parameters
# ─────────────────────────────────────────────
BOX    = 80.0
R_OBS  = 2.0

# pore-size sweep: (n_obstacles, label, color)
PORE_CONFIGS = [
    (120, "Silt (dense)",  PALETTE["silt"]),
    (70,  "Loam (medium)", PALETTE["loam"]),
    (30,  "Sand (sparse)", PALETTE["sand"]),
]

PE_VALUES = [5, 20, 60] if not FAST else [20, 60]
DT        = 1e-3
D_R       = 1.0
RECORD    = 20
N_STEPS   = (200_000 if not FAST else 20_000)
N_SEEDS   = (3 if not FAST else 1)

# polymer params
N_MONO  = 20
KF_FLEX = 30.0
KF_INEX = 500.0
KAPPA_FLEX = 0.0
KAPPA_SEMI = 200.0
N_POLY_STEPS   = (100_000 if not FAST else 15_000)
POLY_RECORD    = 50

print("=" * 55)
print("  Active Transport in Porous Media — Sneha Chakrabarti")
print("=" * 55)
print(f"  Mode: {'FAST' if FAST else 'FULL'}")
print()

# ─────────────────────────────────────────────────────────────
# Helper: build media, run ABP ensemble, return observables
# ─────────────────────────────────────────────────────────────

def run_abp_ensemble(n_obs, Pe, n_seeds=N_SEEDS):
    media = PorousMedia(box_size=BOX, n_obstacles=n_obs, r_obs=R_OBS, seed=7)
    all_pos, all_labels = [], []
    for s in range(n_seeds):
        sim = ABPSimulation(Pe=Pe, D_R=D_R, media=media, dt=DT, seed=s*13+1)
        pos, _, _ = sim.run(N_STEPS, record_every=RECORD)
        labels = classify_frames(pos, window=20, threshold=0.6,
                                 v_thresh=0.03, dt_rec=DT*RECORD)
        all_pos.append(pos)
        all_labels.append(labels)
    # ensemble MSD
    dt_rec = DT * RECORD
    lags, msd0 = compute_msd(all_pos[0], max_lag_fraction=0.3)
    msd_ens = msd0.copy()
    for pos in all_pos[1:]:
        _, m = compute_msd(pos, max_lag_fraction=0.3)
        msd_ens += m
    msd_ens /= len(all_pos)
    p_trap = np.mean([trapping_probability(l) for l in all_labels])
    D_eff  = effective_diffusivity(lags, msd_ens, dt_rec)
    lc, _  = crossover_length(lags, msd_ens, dt_rec)
    return lags * dt_rec, msd_ens, p_trap, D_eff, lc, media

# ─────────────────────────────────────────────────────────────
# FIG 1 — Porous media visualisation + chord-length distribution
# ─────────────────────────────────────────────────────────────
print("→ Fig 1: Porous media + chord-length distributions")

fig, axes = plt.subplots(1, 4, figsize=(14, 3.5))
chord_data = {}
for idx, (n_obs, label, color) in enumerate(PORE_CONFIGS):
    media = PorousMedia(box_size=BOX, n_obstacles=n_obs, r_obs=R_OBS, seed=7)
    ax = axes[idx]
    # draw box
    for cx, cy in media.centers:
        circle = plt.Circle((cx, cy), R_OBS, color=color, alpha=0.55)
        ax.add_patch(circle)
    ax.set_xlim(0, BOX); ax.set_ylim(0, BOX)
    ax.set_aspect('equal')
    ax.set_title(label, fontsize=10, color=color)
    ax.axis('off')
    L = media.characteristic_pore_length(n_chords=(2000 if not FAST else 500))
    chord_data[label] = (color, L)
    ax.text(0.05, 0.05, f"L ≈ {L:.1f} σ", transform=ax.transAxes,
            fontsize=9, color=color)

# chord-length exponential fits
ax = axes[3]
x  = np.linspace(0, 35, 200)
for label, (color, L) in chord_data.items():
    ax.plot(x, expon.pdf(x, scale=L), color=color, label=f"{label}\nL={L:.1f}σ")
ax.set_xlabel("Chord length ℓ (σ)")
ax.set_ylabel("P(ℓ)")
ax.set_title("Pore-size distributions", fontsize=10)
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f"{RESULTS}/fig1_porous_media.png", bbox_inches='tight')
plt.close()
print("   saved fig1_porous_media.png")

# ─────────────────────────────────────────────────────────────
# FIG 2 — MSD for 3 Pe × 3 pore configs
# ─────────────────────────────────────────────────────────────
print("→ Fig 2: MSD sweep Pe × pore size")

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=False)
pe_colors  = ["#1A535C", "#FF6B35", "#9B2335"]
linestyles = ["-", "--", ":"]

msd_store = {}   # (n_obs, Pe) → (t, msd)

for col_idx, (n_obs, label, pcolor) in enumerate(PORE_CONFIGS):
    ax = axes[col_idx]
    for pi, Pe in enumerate(PE_VALUES):
        print(f"     MSD: {label}, Pe={Pe}")
        t, msd, _, _, _, _ = run_abp_ensemble(n_obs, Pe, n_seeds=N_SEEDS)
        msd_store[(n_obs, Pe)] = (t, msd)
        ax.loglog(t, msd, color=pe_colors[pi % 3],
                  ls=linestyles[pi % 3], label=f"Pe={Pe}")
    # guide lines
    t_ref = np.array([t[5], t[len(t)//2], t[-1]])
    ax.loglog(t_ref, 0.5*t_ref**2, 'k:', alpha=0.3, lw=1, label='~τ²')
    ax.loglog(t_ref, 0.8*t_ref,    'k--', alpha=0.3, lw=1, label='~τ')
    ax.set_title(label, color=pcolor)
    ax.set_xlabel("Lag time τ (τ)")
    if col_idx == 0:
        ax.set_ylabel("MSD (σ²)")
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f"{RESULTS}/fig2_msd_sweep.png", bbox_inches='tight')
plt.close()
print("   saved fig2_msd_sweep.png")

# ─────────────────────────────────────────────────────────────
# FIG 3 — D_eff vs L  (for each Pe)
# ─────────────────────────────────────────────────────────────
print("→ Fig 3: Deff vs L")

fig, ax = plt.subplots(figsize=(6, 4))
n_obs_list = [n for n, _, _ in PORE_CONFIGS]
L_list = []
for n_obs in n_obs_list:
    m = PorousMedia(box_size=BOX, n_obstacles=n_obs, r_obs=R_OBS, seed=7)
    L_list.append(m.characteristic_pore_length(n_chords=(1500 if not FAST else 400)))

for pi, Pe in enumerate(PE_VALUES):
    Deff_list = []
    for n_obs in n_obs_list:
        t, msd = msd_store[(n_obs, Pe)]
        Deff_list.append(effective_diffusivity(np.arange(len(msd)), msd, t[1]-t[0]))
    ax.plot(L_list, Deff_list, 'o-', color=pe_colors[pi % 3],
            label=f"Pe={Pe}", ms=8)

ax.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.6, label='Bulk limit')
ax.set_xlabel("Characteristic pore length L (σ)")
ax.set_ylabel(r"$D_{\rm eff}$ (σ²/τ)")
ax.set_title("Effective diffusivity vs pore size")
ax.legend()
plt.tight_layout()
plt.savefig(f"{RESULTS}/fig3_deff_vs_L.png", bbox_inches='tight')
plt.close()
print("   saved fig3_deff_vs_L.png")

# ─────────────────────────────────────────────────────────────
# FIG 4 — P_trapped vs L + inverse spreading rate (Datta-style)
# ─────────────────────────────────────────────────────────────
print("→ Fig 4: Trapping probability vs L")

fig, ax = plt.subplots(figsize=(6, 4))
ax2 = ax.twinx()

Pe_trap = PE_VALUES[-1]   # use highest Pe
ptrap_list, deff_inv_list = [], []
for n_obs in n_obs_list:
    t, msd, ptrap, Deff, lc, _ = run_abp_ensemble(n_obs, Pe_trap, n_seeds=N_SEEDS)
    ptrap_list.append(ptrap)
    deff_inv_list.append(1.0 / max(Deff, 1e-4))

ax.plot(L_list, ptrap_list, 's--', color=PALETTE["silt"],
        ms=9, label=r"$P_{\rm trapped}$")
ax2.plot(L_list, deff_inv_list, '*-', color=PALETTE["sand"],
         ms=11, label=r"$1/D_{\rm eff}$")

ax.set_xlabel("Characteristic pore length L (σ)")
ax.set_ylabel(r"$P_{\rm trapped}$", color=PALETTE["silt"])
ax2.set_ylabel(r"$1/D_{\rm eff}$ (τ/σ²)", color=PALETTE["sand"])
ax.set_title(f"Pore confinement fingerprint  (Pe={Pe_trap})")
lines1, lb1 = ax.get_legend_handles_labels()
lines2, lb2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, lb1+lb2, loc='upper right')
plt.tight_layout()
plt.savefig(f"{RESULTS}/fig4_trapping_vs_L.png", bbox_inches='tight')
plt.close()
print("   saved fig4_trapping_vs_L.png")

# ─────────────────────────────────────────────────────────────
# FIG 5 — Crossover length λ_c vs L
# ─────────────────────────────────────────────────────────────
print("→ Fig 5: Crossover length λ_c vs L")

fig, ax = plt.subplots(figsize=(6, 4))
for pi, Pe in enumerate(PE_VALUES):
    lc_list = []
    for n_obs in n_obs_list:
        t, msd = msd_store[(n_obs, Pe)]
        lags_idx = np.arange(len(msd))
        lc, _ = crossover_length(lags_idx, msd, t[1]-t[0])
        lc_list.append(lc)
    ax.plot(L_list, lc_list, 'D-', color=pe_colors[pi % 3],
            ms=8, label=f"Pe={Pe}")

# bulk value
ax.axhline(5.0, color='k', ls='--', lw=1, alpha=0.5, label='Bulk limit')
ax.plot(L_list, L_list, 'k:', lw=1, alpha=0.4, label='λ_c = L')
ax.set_xlabel("Characteristic pore length L (σ)")
ax.set_ylabel(r"Crossover length $\lambda_c$ (σ)")
ax.set_title("Ballistic→diffusive crossover vs pore size")
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(f"{RESULTS}/fig5_crossover_vs_L.png", bbox_inches='tight')
plt.close()
print("   saved fig5_crossover_vs_L.png")

# ─────────────────────────────────────────────────────────────
# FIG 6 — Polymer COM MSD: linear vs ring, flexible vs semiflexible
# ─────────────────────────────────────────────────────────────
print("→ Fig 6: Polymer MSD — topology × stiffness")
media_poly = PorousMedia(box_size=BOX, n_obstacles=60, r_obs=R_OBS, seed=7)

configs_poly = [
    ("Flexible Linear",   "linear", KF_FLEX, KAPPA_FLEX, PALETTE["linear"],   "-"),
    ("Flexible Ring",     "ring",   KF_FLEX, KAPPA_FLEX, PALETTE["ring"],     "--"),
    ("Semiflexible Linear","linear",KF_FLEX, KAPPA_SEMI, PALETTE["flex"],     "-."),
    ("Semiflexible Ring", "ring",   KF_FLEX, KAPPA_SEMI, PALETTE["semi"],     ":"),
]

fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
Pe_poly_list = [10, 60] if not FAST else [60]

for ax_idx, Pe_p in enumerate(Pe_poly_list):
    ax = axes[ax_idx] if len(Pe_poly_list) > 1 else axes[0]
    ax.set_title(f"Pe = {Pe_p}")
    for cfg_label, top, kf, kap, col, ls in configs_poly:
        print(f"     Polymer MSD: {cfg_label}, Pe={Pe_p}")
        poly = ActivePolymer(N=N_MONO, topology=top, Pe=Pe_p, kf=kf,
                              kappa=kap, D_R=D_R, media=media_poly,
                              dt=5e-4, seed=1)
        coms, _ = poly.run(N_POLY_STEPS, record_every=POLY_RECORD)
        lags, msd = compute_msd(coms, max_lag_fraction=0.25)
        dt_rec = 5e-4 * POLY_RECORD
        t = lags * dt_rec
        ax.loglog(t, msd, color=col, ls=ls, label=cfg_label)
    ax.loglog(t[3:], 0.3*t[3:]**2, 'k:', lw=0.8, alpha=0.3)
    ax.loglog(t[len(t)//3:], 0.8*t[len(t)//3:], 'k--', lw=0.8, alpha=0.3)
    ax.set_xlabel("Lag time τ")
    if ax_idx == 0:
        ax.set_ylabel("COM MSD (σ²)")
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f"{RESULTS}/fig6_polymer_msd.png", bbox_inches='tight')
plt.close()
print("   saved fig6_polymer_msd.png")

# ─────────────────────────────────────────────────────────────
# FIG 7 — Mean Rg vs Pe: linear vs ring
# ─────────────────────────────────────────────────────────────
print("→ Fig 7: Mean Rg vs Pe")

Pe_rg_sweep = [0, 5, 10, 20, 40, 60] if not FAST else [0, 20, 60]
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=False)
media_rg = PorousMedia(box_size=BOX, n_obstacles=60, r_obs=R_OBS, seed=7)
xi_pore  = media_rg.characteristic_pore_length(n_chords=(500 if not FAST else 200))

for ax_idx, (kf, kap, stiffness_label) in enumerate([
        (KF_FLEX, KAPPA_FLEX, "Flexible"),
        (KF_FLEX, KAPPA_SEMI, "Semiflexible")]):
    ax = axes[ax_idx]
    for top, col, marker in [("linear", PALETTE["linear"], "o"),
                               ("ring",   PALETTE["ring"],   "s")]:
        print(f"     Rg sweep: {stiffness_label} {top}")
        rg_vals = []
        for Pe_p in Pe_rg_sweep:
            poly = ActivePolymer(N=N_MONO, topology=top, Pe=Pe_p, kf=kf,
                                  kappa=kap, D_R=D_R, media=media_rg,
                                  dt=5e-4, seed=2)
            _, rgs = poly.run(N_POLY_STEPS, record_every=POLY_RECORD)
            rg_vals.append(np.mean(rgs[len(rgs)//3:]) / xi_pore)
        ax.plot(Pe_rg_sweep, rg_vals, f'{marker}-', color=col,
                ms=8, label=f"{top.capitalize()}")
    ax.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.5, label=r'⟨Rg⟩/ξ=1')
    ax.set_xlabel("Péclet number Pe")
    ax.set_ylabel(r"$\langle R_g\rangle / \xi$")
    ax.set_title(f"{stiffness_label} polymers")
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(f"{RESULTS}/fig7_rg_vs_pe.png", bbox_inches='tight')
plt.close()
print("   saved fig7_rg_vs_pe.png")

# ─────────────────────────────────────────────────────────────
# FIG 8 — Asphericity A vs Pe
# ─────────────────────────────────────────────────────────────
print("→ Fig 8: Asphericity vs Pe")

fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
media_asp = PorousMedia(box_size=BOX, n_obstacles=60, r_obs=R_OBS, seed=7)

for ax_idx, (kf, kap, stiffness_label) in enumerate([
        (KF_FLEX, KAPPA_FLEX, "Flexible"),
        (KF_FLEX, KAPPA_SEMI, "Semiflexible")]):
    ax = axes[ax_idx]
    for top, col, marker in [("linear", PALETTE["linear"], "o"),
                               ("ring",   PALETTE["ring"],   "s")]:
        print(f"     Asphericity: {stiffness_label} {top}")
        a_vals = []
        for Pe_p in Pe_rg_sweep:
            poly = ActivePolymer(N=N_MONO, topology=top, Pe=Pe_p, kf=kf,
                                  kappa=kap, D_R=D_R, media=media_asp,
                                  dt=5e-4, seed=3)
            A_list = []
            for s in range(N_POLY_STEPS):
                poly.step()
                if s % POLY_RECORD == 0 and s > N_POLY_STEPS//3:
                    A_list.append(poly.asphericity())
            a_vals.append(np.mean(A_list))
        ax.plot(Pe_rg_sweep, a_vals, f'{marker}-', color=col,
                ms=8, label=f"{top.capitalize()}")
    ax.set_xlabel("Péclet number Pe")
    ax.set_ylabel("Asphericity ⟨A⟩")
    ax.set_title(f"{stiffness_label} — shape anisotropy")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(f"{RESULTS}/fig8_asphericity_vs_pe.png", bbox_inches='tight')
plt.close()
print("   saved fig8_asphericity_vs_pe.png")

# ─────────────────────────────────────────────────────────────
# FIG 9 — Phase diagram: P_trapped(L, Pe)
# ─────────────────────────────────────────────────────────────
print("→ Fig 9: Phase diagram P_trapped(L, Pe)")

Pe_grid   = np.array([5, 10, 20, 40, 60] if not FAST else [10, 40])
n_obs_grid = [120, 90, 70, 50, 30] if not FAST else [120, 70, 30]

L_grid  = []
Pt_grid = np.zeros((len(Pe_grid), len(n_obs_grid)))

for j, n_obs in enumerate(n_obs_grid):
    m  = PorousMedia(box_size=BOX, n_obstacles=n_obs, r_obs=R_OBS, seed=7)
    L_grid.append(m.characteristic_pore_length(n_chords=(500 if not FAST else 200)))
    for i, Pe in enumerate(Pe_grid):
        print(f"     Phase map: Pe={Pe}, n_obs={n_obs}")
        _, _, ptrap, _, _, _ = run_abp_ensemble(n_obs, Pe, n_seeds=1)
        Pt_grid[i, j] = ptrap

fig, ax = plt.subplots(figsize=(7, 4.5))
im = ax.contourf(L_grid, Pe_grid, Pt_grid,
                 levels=20, cmap='RdYlGn_r', vmin=0, vmax=1)
cbar = fig.colorbar(im, ax=ax, label=r"$P_{\rm trapped}$")
ax.set_xlabel("Characteristic pore length L (σ)")
ax.set_ylabel("Péclet number Pe")
ax.set_title("Phase diagram: trapping probability")

# annotate regimes
ax.text(L_grid[0]*1.1, Pe_grid[-1]*0.9, "Trapped\n(slow diffusion)",
        color='white', fontsize=9, ha='left', va='top',
        fontweight='bold')
ax.text(L_grid[-1]*0.85, Pe_grid[0]*1.1, "Run-and-tumble\n(fast diffusion)",
        color='black', fontsize=9, ha='right', va='bottom',
        fontweight='bold')
plt.tight_layout()
plt.savefig(f"{RESULTS}/fig9_phase_diagram.png", bbox_inches='tight')
plt.close()
print("   saved fig9_phase_diagram.png")

# ─────────────────────────────────────────────────────────────
# Summary table
# ─────────────────────────────────────────────────────────────
print()
print("=" * 55)
print("  Summary: Pore length vs. effective diffusivity")
print("  (Pe = {})".format(PE_VALUES[-1]))
print("=" * 55)
print(f"  {'Regime':<18} {'L (σ)':>8} {'D_eff':>12} {'P_trap':>10}")
print("  " + "-"*50)
labels_p = ["Silt (dense)", "Loam (medium)", "Sand (sparse)"]
for n_obs, lbl, _ in PORE_CONFIGS:
    t, msd = msd_store[(n_obs, PE_VALUES[-1])]
    Deff = effective_diffusivity(np.arange(len(msd)), msd, t[1]-t[0])
    _, _, ptrap, _, _, m = run_abp_ensemble(n_obs, PE_VALUES[-1], n_seeds=1)
    L = m.characteristic_pore_length(n_chords=500)
    print(f"  {lbl:<18} {L:>8.2f} {Deff:>12.4f} {ptrap:>10.3f}")

print()
print("  All figures saved in ./results/")
print("=" * 55)

