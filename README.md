# pore-geometry-as-ecological-filter-ABP-and-active-polymer-transport-2D

How does pore-scale geometry act as a physical filter on active transport: quantifying the confinement-driven transition from run-and-tumble to trapped motion in single ABPs across silt-to-sand pore regimes.

---

# Active Transport in Porous Media  
### A Computational Physics Project  

**Author:** Sneha Chakrabarti  
**Affiliation:** IISER Kolkata (3rd year BS-MS, Physics Major)  
**Contact:** sc23ms114@iiserkol.ac.in  

---

## Overview

This project investigates how pore geometry controls active transport at the level of:

- **Single Active Brownian Particles (ABPs)** — confinement shifts motion from run-and-tumble to trapped states, modifying effective diffusivity and trapping probability.

This connects to:
- Bhattacharyay (2011)  
- Al Harraq et al., PRX Life (2026)  

---

## Scientific Background

### Active Brownian Particles (ABPs)

dx/dt = v₀ ê(θ) + f_obs + √(2 k_B T) η(t)  
dθ/dt = √(2 D_R) ξ(t)

Pe = F_a σ / (k_B T)  
τ_R = 1 / D_R  

---

### Pore Size Characterization

P(ℓ) = (1 / L) exp(−ℓ / L)

L = characteristic pore size  

---

### Run Classification

δ(t) = net displacement / path length  

trapped if δ < 0.6 or speed < threshold  

---

### MSD Scaling

MSD(τ) ~ τ² (ballistic)  
MSD(τ) ~ τ (diffusive)  

λ_c = √MSD(τ_c)

---

## Repository Structure

active-transport-porous-media/  
├── src/  
├── results/  
│   ├── fig1_porous_media.png  
│   ├── fig2_msd_sweep.png  
│   ├── fig3_deff_vs_L.png  
│   ├── fig4_trapping_vs_L.png  
│   └── fig5_crossover_vs_L.png  
├── docs/  
└── README.md  

---

## Installation

pip install numpy matplotlib scipy  

Python ≥ 3.9 recommended.

---

## Quick Run

cd src  
python run_analysis.py --fast  

---

## Results Summary

### Fig 1 — Porous Media Characterisation

P(ℓ) = (1 / L) exp(−ℓ / L)

| Texture | N obstacles | L (σ) | Interpretation |
|--------|------------|------|----------------|
| Silt   | 120        | 10.6 | Highly confining |
| Loam   | 70         | 19.0 | Marginal |
| Sand   | 30         | 36.7 | Nearly free |

Inference: exponential pore statistics.

---

### Fig 2 — MSD vs Lag Time

MSD ~ τ² (short times)  
MSD ~ τ (long times)  

ℓ_p = v₀ / D_R = Pe  

Inference:
- MSD suppressed in silt vs sand  
- strongest suppression at low Pe  

---

### Fig 3 — Effective Diffusivity vs L

D_eff = MSD(τ) / (4 τ)

| Texture | Pe | L (σ) | D_eff | Suppression |
|--------|----|------|-------|------------|
| Silt   | 20 | 10.6 | 126.0 | 1.8× |
| Loam   | 20 | 19.0 | 176.1 | 1.3× |
| Sand   | 20 | 36.7 | 226.3 | — |
| Silt   | 60 | 10.6 | 139.0 | 1.4× |
| Loam   | 60 | 19.0 | 174.0 | 1.1× |
| Sand   | 60 | 36.7 | 200.0 | — |

Inference:
- D_eff increases with L  
- confinement persists even at high activity  

---

### Fig 4 — Trapping Probability vs L

δ(t) = Δr_net / Δr_path  

P_trapped = trapped_frames / total_frames  

| Texture | L (σ) | P_trapped | 1/D_eff |
|--------|------|-----------|----------|
| Silt   | 10.6 | 0.113     | 7.94e−3 |
| Loam   | 19.0 | 0.076     | 5.68e−3 |
| Sand   | 36.7 | 0.049     | 4.42e−3 |

Key relation:

D_eff ∝ (1 − P_trapped)

Inference:
- trapping directly controls diffusion  

---

### Fig 5 — Crossover Length λ_c vs L

λ_c = √MSD(τ_c)

| Texture | L (σ) | λ_c (Pe=20) | λ_c (Pe=60) |
|--------|------|-------------|-------------|
| Silt   | 10.6 | 9.9         | 25.0        |
| Loam   | 19.0 | 10.0        | 33.3        |
| Sand   | 36.7 | 8.0         | 27.8        |

Regimes:

Boundary-limited (L < ℓ_p):  
λ_c ∝ L  

Motility-limited (L > ℓ_p):  
λ_c → ℓ_p  

Phase boundary:  
L* = ℓ_p = Pe / D_R  

---

## Physical Interpretation

Particles alternate between running and trapped states.

D_eff ∝ (1 − P_trapped)

Transport is controlled by time spent trapped.

---

## Connection to Soil Ecology

- Small pores → trapping dominates  
- Large pores → free swimming  
- Explains failure of chemotaxis in silt  

---

## Limitations

- 2D simulations  
- No hydrodynamics  
- No chemotaxis  
- Monodisperse obstacles  

---

## Future Directions

- Extend to 3D  
- Add chemotaxis  
- Multi-particle effects  

---

## References

- Bhattacharyay (2011)  
- Al Harraq et al., PRX Life (2026)  
- Torquato & Lu (1993)  

---

## Citation

Chakrabarti, S. (2025). Active Transport in Porous Media: ABP Simulations. GitHub repository. IISER Kolkata.
