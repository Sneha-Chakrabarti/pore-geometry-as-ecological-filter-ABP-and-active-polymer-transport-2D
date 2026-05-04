# Active Transport in Porous Media
### A Computational Physics Project

**Author:** Sneha Chakrabarti  
**Affiliation:** IISER Kolkata (3rd year BS-MS, Physics Major)  
**Contact:** sc23ms114@iiserkol.ac.in

---

## Overview

This project investigates how **pore geometry controls active transport** at two levels:

1. **Single Active Brownian Particles (ABPs)** — how confinement shifts run-and-tumble dynamics to trapped motion, modulating effective diffusivity and trapping probability as a function of pore size and Péclet number.

2. **Active polymer chains and rings** — how topology (linear vs. ring) and stiffness (flexible vs. semiflexible) interact with pore confinement to determine conformational swelling/shrinking and COM diffusion.

This directly bridges three recent papers:
- Bhattacharyay (2011) — directed transport and symmetry-breaking in equilibrium
- Theeyancheri, Chaki, Bhattacharjee, Chakrabarti, JCP 159, 014902 (2023) — active polymer chains/rings in porous media
- Al Harraq, Choi, Datta, Shaevitz, PRX Life 4, 013034 (2026) — bacterial motility & chemotaxis regulated by soil texture

---

## Scientific Background

### Active Brownian Particles (ABPs)

An ABP is a minimal model for a self-propelled microswimmer (bacterium, Janus particle). In overdamped Langevin dynamics:

```
dx/dt = v₀ ê(θ) + f_obs + √(2kBT) η(t)
dθ/dt = √(2D_R) ξ(t)
```

where:
- `v₀ = Pe × kBT/(γσ)` — self-propulsion speed (Pe = Péclet number)  
- `ê(θ)` — unit orientation vector  
- `D_R` — rotational diffusion coefficient; persistence time `τ_R = 1/D_R`  
- `f_obs` — WCA repulsion from obstacles  
- `η, ξ` — Gaussian white noise

**Key dimensionless number:** Pe = Fa·σ/(kBT) — ratio of active to thermal forces.

### Pore-Size Characterization

Pore geometry is quantified via the **chord-length distribution** (Torquato & Lu 1993):

```
P(ℓ) = (1/L) exp(−ℓ/L)
```

The characteristic length `L` is the mean pore size estimated from random line-chord sampling. This mirrors the experimental method of Al Harraq et al. 2026 (cryolite soil microcosms) and maps to:
- Silt: small `L` → highly confining
- Loam: intermediate `L`
- Sand: large `L` → weakly confining

### Run Classification

Each trajectory frame is classified as **running** or **trapped** using a displacement-ratio algorithm:

```
δ(t) = net_displacement / path_length   (over window w)
trapped if δ < 0.6  OR  instantaneous_speed < v_thresh
```

This reproduces the classification of Al Harraq et al. PRX Life 2026.

### MSD Scaling Regimes

```
MSD(τ) ~ τ²    (ballistic, short lag time — persistent swimming)
MSD(τ) ~ τ     (diffusive, long lag time — tumbling + pore collisions)
```

The crossover length `λ_c = √MSD(τ_c)` is extracted by segmented log-log regression. In small pores, `λ_c ∝ L` (boundary-limited); for large pores, `λ_c` saturates to the bulk value (motility-limited).

### Active Polymer Model

The polymer is a chain of N monomers with:
- **FENE bond** (Kremer-Grest): `V_FENE = -k_f r²_max/2 × ln[1 − (r/r_max)²]`
- **WCA excluded volume**: `V_WCA = 4ε[(σ/r)¹² − (σ/r)⁶] + ε`
- **Bending rigidity** (Kratky-Porod): `V_bend = κ(1 − cos φ_i)`
- **Active force** on each monomer: `Fa ê(θᵢ)` with independent rotational diffusion

**Topology:**
- Linear: N monomers, N−1 bonds — can adopt rod-like conformations, navigates pores readily
- Ring: N monomers, N bonds (terminal connected) — closed loop geometry causes trapping in confining pores

**Shape descriptor** — Asphericity from gyration tensor eigenvalues λ₁ ≤ λ₂:
```
A = (λ₂ − λ₁)² / (λ₁ + λ₂)²
A = 0: circular    A = 1: fully extended rod
```

---

## Repository Structure

```
active-transport-porous-media/
├── src/
│   ├── abp_porous.py        # Core ABP simulation + analysis utilities
│   ├── active_polymer.py    # Active polymer (linear/ring) simulation
│   └── run_analysis.py      # Master script: runs all simulations, saves figures
├── results/                 # Auto-generated figures (created on run)
│   ├── fig1_porous_media.png
│   ├── fig2_msd_sweep.png
│   ├── fig3_deff_vs_L.png
│   ├── fig4_trapping_vs_L.png
│   ├── fig5_crossover_vs_L.png
│   ├── fig6_polymer_msd.png
│   ├── fig7_rg_vs_pe.png
│   ├── fig8_asphericity_vs_pe.png
│   └── fig9_phase_diagram.png
├── docs/
│   └── theory_notes.md      # Detailed derivations and physical interpretation
├── requirements.txt
└── README.md
```

---

## Installation & Usage

### Requirements

```bash
pip install numpy matplotlib scipy
```
Python ≥ 3.9 recommended.

### Quick Run (reduced statistics, ~3–5 min)

```bash
cd src
python run_analysis.py --fast
```

### Full Run (production quality, ~30–60 min)

```bash
cd src
python run_analysis.py
```

All figures are saved to `results/`.

### Run Individual Modules

```python
from abp_porous import PorousMedia, ABPSimulation, compute_msd

# Build porous medium
media = PorousMedia(box_size=80.0, n_obstacles=70, r_obs=2.0, seed=42)
L = media.characteristic_pore_length()
print(f"Characteristic pore length L = {L:.2f} σ")

# Run ABP simulation
sim = ABPSimulation(Pe=20, D_R=1.0, media=media, dt=1e-3, seed=0)
positions, wrapped, thetas = sim.run(n_steps=100_000, record_every=20)

# Compute MSD
lags, msd = compute_msd(positions)
```

```python
from active_polymer import ActivePolymer

# Semiflexible ring in porous medium
ring = ActivePolymer(
    N=20, topology='ring',
    Pe=60, kf=30.0, kappa=200.0,
    D_R=1.0, media=media, dt=5e-4
)
com_positions, rg_trajectory = ring.run(n_steps=200_000, record_every=50)
print(f"Mean Rg = {rg_trajectory.mean():.2f} σ")
print(f"Mean asphericity = {ring.asphericity():.3f}")
```

---

## Results Summary

### Fig 1 — Porous Media Visualization
Three packing densities (silt/loam/sand analog) with exponential chord-length distributions. The characteristic length L spans from ~5σ (dense/silt) to ~20σ (sparse/sand).

### Fig 2 — MSD vs Lag Time
- **Ballistic regime** (τ²) at short lag times: active swimming dominates
- **Diffusive regime** (τ¹) at long lag times: pore collisions randomize direction
- Higher Pe → larger MSD amplitude but same qualitative scaling
- Silt (dense) suppresses MSD amplitude by ~10× relative to sand (sparse)

### Fig 3 — Effective Diffusivity vs L
D_eff increases monotonically with L. In silt-like confinement, D_eff drops by roughly an order of magnitude relative to the bulk value. At large Pe, D_eff is enhanced but the suppression in tight pores persists — **confinement dominates over activity at small L**.

### Fig 4 — Trapping Probability vs L
P_trapped decreases monotonically with increasing pore size and inversely mirrors 1/D_eff, confirming that **microscopic trapping events directly set the macroscopic transport rate** (Al Harraq et al. 2026).

### Fig 5 — Crossover Length λ_c vs L
- For L < 20σ: λ_c ∝ L (boundary-limited regime — particles hit walls before tumbling)
- For L > 20σ: λ_c saturates to bulk value (motility-limited regime — intrinsic tumbling governs)
This crossover matches the experimental observation in cryolite soils (Al Harraq et al. 2026, Fig. 3e).

### Fig 6 — Polymer MSD: Topology × Stiffness
- Flexible linear and ring chains both diffuse, with linear chains migrating faster
- Semiflexible linear chains navigate smoothly via rod-like conformations
- **Semiflexible rings exhibit extended subdiffusion** — circular topology causes trapping in pore confinements (Theeyancheri et al. 2023)

### Fig 7 — Mean Rg vs Pe
- **Flexible chains**: activity-induced swelling for both linear and ring
- **Semiflexible rings**: activity-induced shrinking — pore collisions push monomers inward (ring cannot extend rod-like)
- **Semiflexible linear**: shrinks at low Pe, swells at high Pe (consistent with Theeyancheri et al. 2023, Fig. 6)

### Fig 8 — Asphericity A vs Pe
Asphericity increases with activity for all architectures. Semiflexible linear chains have higher A (more rod-like) than flexible ones. Semiflexible rings show lower A, consistent with their circular geometry and trapping.

### Fig 9 — Phase Diagram P_trapped(L, Pe)
The 2D phase space reveals:
- **Trapped phase** (upper left): small L, high Pe — particles push hard into walls and get stuck
- **Run-and-tumble phase** (lower right): large L, low Pe — free swimming dominates
- The boundary between phases shifts to smaller L as Pe increases

---

## Physical Interpretation

### Why rings trap more than linear chains

A semiflexible ring's closed-loop geometry forces it to adopt a circular conformation. When this circle encounters a pore throat, no part of the chain can snake through independently — the entire ring must deform significantly, which costs bending energy ∝ κ. By contrast, a linear chain can insert one end like a needle into the pore. This topology-controlled trapping is the physical mechanism behind Fig. 6 and 7, and explains why ring-topology drug carriers have longer blood circulation times (Chen et al. 2009).

### Why P_trapped ∝ 1/D_eff

The trapping probability P_trapped controls how long a particle is immobilized per unit time. The effective diffusivity of a particle that alternates between running (speed v₀, run length λ) and trapped phases (duration τ_trap) scales as:

```
D_eff ~ (v₀ λ) × (1 − P_trapped) / (1 + P_trapped × τ_trap / τ_run)
```

So D_eff is a monotonically decreasing function of P_trapped — the correspondence in Fig. 4 is a direct consequence of this relation.

### Connection to soil ecology (Datta lab 2026)

In the experimental system of Al Harraq et al. 2026, the characteristic pore length L determines whether E. coli can reach Arabidopsis roots via chemotaxis. Our simulations quantitatively reproduce:
- The order-of-magnitude reduction in D_eff from sandy (L~48μm) to silty (L~7μm) soils
- The monotonic increase of P_trapped as L decreases
- The saturation of λ_c at bulk values for large L

This suggests that the **physical confinement mechanism** (not chemical gradient suppression) dominates chemotactic failure in silty soils.

---

## Limitations and Future Directions

| Limitation | Future Extension |
|---|---|
| 2D simulations only | 3D porous media (Theeyancheri et al. mention qualitative changes) |
| Overdamped dynamics | Full inertial Langevin (relevant for Bhattacharyay 2011 dimer model) |
| Monodisperse obstacles | Polydisperse pore sizes (realistic soil) |
| No hydrodynamic interactions | Include Oseen tensor for dense suspensions |
| Single-particle | Multi-particle: collective effects, jamming, MIPS |
| Passive obstacles | Active/deformable obstacles (biological tissues) |
| No chemical gradients | Chemotaxis: ABP + chemical field (Bhattacharjee & Datta 2021) |
| Symmetric activity | Asymmetric damping → directed transport (Bhattacharyay 2011) |

The last point is particularly motivated by Bhattacharyay (2011): a symmetry-broken dimer with **different damping constants** on each particle can undergo steady directed transport in equilibrium. An active version of this — where activity replaces thermal motion — is an open problem directly relevant to understanding intracellular motor proteins.

---

## Relation to Ongoing Research Interests

This project is designed to strengthen research in:

1. **Statistical mechanics of active matter** — connects to Dr. Rumi De's lab (predator-prey, cooperative dynamics)
2. **Defect dynamics in active media** — the trapping/escaping phenomenology of semiflexible rings parallels defect-mediated dynamics in active nematics (upcoming IISER Pune project)
3. **Active transport in porous environments** — directly relevant to planned IIT Bombay project with Prof. Rajarshi Chakrabarti
4. **Soft condensed matter** — Rg scaling, polymer conformation analysis align with planned 4th-year project with Prof. Arindam Kundagrami

---

## References

1. A. Bhattacharyay, "Directed transport in equilibrium," arXiv:1008.4992 (2011)
2. L. Theeyancheri, S. Chaki, T. Bhattacharjee, R. Chakrabarti, "Active dynamics of linear chains and rings in porous media," *J. Chem. Phys.* **159**, 014902 (2023)
3. A. Al Harraq, G. Choi, S. S. Datta, J. W. Shaevitz, "Soil texture regulates bacterial motility and chemotactic recruitment to plant roots," *PRX Life* **4**, 013034 (2026)
4. T. Bhattacharjee and S. S. Datta, "Bacterial hopping and trapping in porous media," *Nat. Commun.* **10**, 2075 (2019)
5. S. Torquato and B. Lu, "Chord-length distribution function for two-phase random media," *Phys. Rev. E* **47**, 2950 (1993)
6. K. Kremer and G. S. Grest, "Dynamics of entangled linear polymer melts," *J. Chem. Phys.* **92**, 5057 (1990)
7. C. Bechinger et al., "Active particles in complex and crowded environments," *Rev. Mod. Phys.* **88**, 045006 (2016)

---

## Citation

If you use this code, please cite:
```
Chakrabarti, S. (2025). Active Transport in Porous Media: 
ABP and Active Polymer Simulations. GitHub repository.
IISER Kolkata.
```
