# pore-geometry-as-ecological-filter-ABP-and-active-polymer-transport-2D
# Active Transport in Porous Media  
### A Computational Physics Project  

**Author:** Sneha Chakrabarti  
**Affiliation:** IISER Kolkata (3rd year BS-MS, Physics Major)  
**Contact:** sc23ms114@iiserkol.ac.in  

This is an independent computational physics project investigating how pore-scale geometry acts as a physical filter on active transport.

The work combines:
- stochastic dynamics (Langevin equations)
- statistical physics of transport
- geometrical characterization of porous media  

and connects directly to recent experimental work on bacterial motility in soil environments.
How does pore-scale geometry act as a physical filter on active transport?

This project quantifies the transition from run-and-tumble to trapped motion in active Brownian particles (ABPs) across silt-to-sand pore regimes.

---

## ⭐ Key Result

Active transport in porous media is governed by a single dimensionless parameter:

L / ℓ_p

- L = characteristic pore size  
- ℓ_p = persistence length = v₀ / D_R  

A transition at L ~ ℓ_p separates:

- Boundary-limited regime → trapping-dominated transport  
- Motility-limited regime → free diffusion  

The effective diffusivity satisfies:

D_eff ∝ (1 − P_trapped)

linking microscopic trapping dynamics directly to macroscopic transport.

---


## 🔬 Model

Active Brownian Particle dynamics:

dx/dt = v₀ ê(θ) + f_obs + √(2 k_B T) η(t)  
dθ/dt = √(2 D_R) ξ(t)

Pe = F_a σ / (k_B T)  
τ_R = 1 / D_R  

---

## 🧠 Key Ingredients

### Pore Geometry

P(ℓ) = (1 / L) exp(−ℓ / L)

L = characteristic pore size  

---

### Run Classification

δ(t) = net displacement / path length  

Trapped if:
- δ < 0.6  
- or instantaneous speed < threshold  

---

### Transport Observables

MSD(τ) ~ τ² (ballistic)  
MSD(τ) ~ τ (diffusive)  

λ_c = √MSD(τ_c)  

D_eff = MSD(τ) / (4 τ)

---

## 📊 Key Results

### MSD suppression across pore sizes
![MSD](results/fig2_msd_sweep.png)

---

### Trapping controls transport
![Trapping](results/fig4_trapping_vs_L.png)

---

### Phase diagram of transport regimes
![Phase](results/fig9_phase_diagram.png)

---

## Results Summary

### Fig 1 — Porous Media Characterisation

P(ℓ) = (1 / L) exp(−ℓ / L)

| Texture | N obstacles | L (σ) | Interpretation |
|--------|------------|------|----------------|
| Silt   | 120        | 10.6 | Highly confining |
| Loam   | 70         | 19.0 | Marginal |
| Sand   | 30         | 36.7 | Nearly free |

---

### Fig 2 — MSD vs Lag Time

- MSD ~ τ² (short times)  
- MSD ~ τ (long times)  

ℓ_p = v₀ / D_R = Pe  

---

### Fig 3 — Effective Diffusivity vs L

D_eff = MSD(τ) / (4 τ)

- D_eff increases with L  
- Confinement persists even at high activity  

---

### Fig 4 — Trapping Probability vs L

δ(t) = Δr_net / Δr_path  

P_trapped = trapped_frames / total_frames  

D_eff ∝ (1 − P_trapped)

---

### Fig 5 — Crossover Length λ_c vs L

λ_c = √MSD(τ_c)

Boundary-limited (L < ℓ_p):  
λ_c ∝ L  

Motility-limited (L > ℓ_p):  
λ_c → ℓ_p  

Phase boundary:  
L* = ℓ_p = Pe / D_R  

---

## 🌍 Why this matters

This model explains how soil pore geometry regulates microbial transport:

- Small pores → trapping dominates → suppressed diffusion  
- Large pores → persistent motion → effective migration  

Provides a physical explanation for experimentally observed limits of bacterial transport in dense soils (Datta Lab, PRX Life 2026).

---

## ⚙️ Repository Structure

active-transport-porous-media/  
├── src/  
├── results/  
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

## Limitations

- 2D simulations  
- No hydrodynamics  
- No chemotaxis  
- Monodisperse obstacles  

---

## Future Directions

- Extend to 3D  
- Add chemotaxis  
- Multi-particle interactions  

---

## References

- Bhattacharyay (2011)  
- Al Harraq et al., PRX Life (2026)  
- Torquato & Lu (1993)  

---

## Citation

Chakrabarti, S. (2025). Active Transport in Porous Media: ABP Simulations. GitHub repository. IISER Kolkata.
