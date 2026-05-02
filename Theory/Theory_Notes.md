
Theory Notes: Active Transport in Porous Media
Sneha Chakrabarti | IISER Kolkata
---
1. Active Brownian Particles — The Overdamped Langevin Picture
1.1 Equations of Motion
For a spherical microswimmer of diameter σ, mass m, in a viscous solvent at temperature T, the full Langevin equation is:
```
m ẍ = -γ ẋ + Fa ê(θ) + f_obs(x) + √(2γkBT) η(t)
```
In the overdamped limit (low Reynolds number, inertia negligible, valid for bacteria):
```
γ ẋ = Fa ê(θ) + f_obs(x) + √(2γkBT) η(t)
θ̇  = √(2DR) ξ(t)
```
Setting γ = 1, kBT = 1, σ = 1 (reduced units):
```
ẋ = v₀ ê(θ) + f_obs + √2 η(t)
θ̇ = √(2DR) ξ(t)
```
where `⟨η_i(t)η_j(t')⟩ = δᵢⱼδ(t−t')`, `⟨ξ(t)ξ(t')⟩ = δ(t−t')`.
1.2 The Péclet Number
The Péclet number Pe = Fa σ / (kBT) = v₀ / (kBT/γσ) measures the ratio of active advection to thermal diffusion. In reduced units:
```
Pe = v₀ = Fa    (since γ = kBT = σ = 1)
```
Pe ≪ 1 : thermally dominated, nearly passive Brownian motion
Pe ~ 1  : active and thermal forces comparable
Pe ≫ 1  : strongly active, long persistent runs
For E. coli: v₀ ≈ 20–30 μm/s, D_R ≈ 0.06 rad²/s, giving Pe ~ 10–30 (motility buffer).
1.3 MSD in Free Space
In free 2D space (no obstacles), the MSD of a single ABP is exactly:
```
MSD(τ) = 4kBT/γ · τ  +  v₀² τ_R² [2τ/τ_R + 2e^{-τ/τ_R} − 2]
```
where τ_R = 1/D_R is the persistence time.
Short times (τ ≪ τ_R): MSD ≈ v₀² τ²  (ballistic, slope 2 on log-log)
Long times (τ ≫ τ_R): MSD ≈ 4D_eff τ  (diffusive, slope 1 on log-log)
with D_eff = kBT/γ + v₀² τ_R / 2 = D₀ + v₀²/(2D_R)
The crossover occurs at τ ~ τ_R, so λ_c ~ v₀/D_R (the persistence length ℓ_p = v₀/D_R).
1.4 Effect of Confinement (Porous Media)
When characteristic pore size L < ℓ_p = v₀/D_R, the particle reaches a wall before it tumbles. The effective run is truncated:
```
λ_c(L) ≈ min(ℓ_p, L)
```
L < ℓ_p: boundary-limited regime, λ_c ∝ L (linear scaling)
L > ℓ_p: motility-limited regime, λ_c → ℓ_p (plateau)
This is precisely what Fig. 5 shows. The crossover between regimes occurs at L* ≈ ℓ_p = v₀/D_R. In our reduced units with D_R = 1: L* ≈ v₀ = Pe.
1.5 Effective Diffusivity in Porous Media
A simple two-state model (run + trap) gives:
```
D_eff = D_run × (1 − P_trapped) / (1 + P_trapped × τ_trap/τ_run)
```
where D_run ~ v₀² τ_R/2. This is why Fig. 3 and Fig. 4 are mirror images of each other.
More precisely, using continuous-time random walk theory:
```
D_eff = λ_c² / (2τ_c) × (1 − P_trapped)
```
where τ_c is the crossover time. Confinement simultaneously reduces λ_c AND increases P_trapped, so D_eff is doubly suppressed in tight pores.
---
2. Pore-Space Geometry and Chord-Length Statistics
2.1 Random Chord Sampling
For a random two-phase medium (solid + void), the chord-length distribution P(ℓ) describes the distribution of straight-line segment lengths through the void phase.
For a Poisson random medium (randomly placed, non-overlapping obstacles):
```
P(ℓ) = (1/L) exp(−ℓ/L)
```
where L is the mean chord length (characteristic pore size). The exponential form arises from the memoryless statistics of random granular packings (Torquato & Lu 1993).
2.2 Maximum-Likelihood Estimate of L
For N observed chord lengths {ℓ₁, ..., ℓ_N}, the MLE for an exponential distribution is:
```
L̂ = (1/N) Σᵢ ℓᵢ = ⟨ℓ⟩
```
This is what `PorousMedia.characteristic_pore_length()` computes.
2.3 Mapping to Soil Textures
The experimental cryolite soils (Al Harraq et al. 2026) have:
Silt:  L ≈ 7–10 μm
Loam:  L ≈ 20–28 μm
Sand:  L ≈ 42–48 μm
E. coli run length ℓ_p ≈ 20–30 μm. So:
Silt: L < ℓ_p → boundary-limited → frequent trapping
Sand: L > ℓ_p → motility-limited → free run-and-tumble
Our simulation uses σ as the unit; with N_mono = 60–120 obstacles in L=80σ, we achieve L ~ 5–15σ, which maps to the same physical regime when σ ~ E. coli body length (~2 μm).
---
3. WCA Potential — Excluded Volume Interaction
The Weeks-Chandler-Andersen (WCA) potential is a purely repulsive, short-range potential:
```
V_WCA(r) = 4ε [(σ_eff/r)¹² − (σ_eff/r)⁶] + ε    if r < 2^{1/6} σ_eff
           = 0                                         if r ≥ 2^{1/6} σ_eff
```
where σ_eff = (σ_particle + σ_obstacle)/2 for particle-obstacle interactions.
Force: `F = −∇V = 4ε/r² [12(σ/r)¹² − 6(σ/r)⁶] r̂`
This acts as a soft elastic wall — particles cannot penetrate obstacles but don't stick to them either (no attractive well). This is appropriate for sterically-stabilized particles and bacteria (hydrophobic BSA coating in the experiments).
---
4. Active Polymer Model
4.1 FENE Bond Potential
The Finitely Extensible Nonlinear Elastic (FENE) potential models covalent bonds with a maximum extension r_max:
```
V_FENE(r) = −(k_f r²_max / 2) × ln[1 − (r/r_max)²]    if r < r_max
           = ∞                                            if r ≥ r_max
```
Force: `F = −k_f r / [1 − (r/r_max)²]`
This diverges as r → r_max, preventing bond breaking. In reduced units with r_max = 1.5σ, the equilibrium bond length ≈ 0.96σ.
4.2 Bending Rigidity — Kratky-Porod Model
For a semiflexible polymer, bending is penalized:
```
V_bend = κ Σᵢ (1 − cos φᵢ)
```
where φᵢ is the angle between consecutive bond vectors (bond i and bond i+1).
The persistence length in 2D for a freely-jointed chain with bending stiffness is:
```
ℓ_p = −b / ln(⟨cos φ⟩)  ≈  2κb / kBT    (for small κ in kBT units)
```
where b is the bond length (~σ in our model). With κ = 200 (kBT units), ℓ_p ~ 200σ ≫ chain contour length N·b = 20σ, giving a rod-like polymer. With κ = 0, we have a freely-jointed chain.
4.3 Activity and the Péclet Number for Polymers
Each monomer has independent orientation θᵢ evolving under:
```
dθᵢ/dt = √(2DR) ξᵢ(t)
```
The active force on monomer i is `Fa ê(θᵢ)`. For N monomers with independent orientations, the net active force on the COM is:
```
F_active,COM = Fa/N × |Σᵢ ê(θᵢ)|
```
At short times (t ≪ τ_R), all orientations are correlated → net force ~ Fa (directed).  
At long times (t ≫ τ_R), orientations decorrelate → net COM force → 0, but fluctuations persist.
4.4 Topology Effect: Linear vs Ring
Linear chain:
N monomers, N−1 FENE bonds
Can elongate along one axis → high asphericity A ≈ 1 for stiff chains
In pores: can insert end-first, navigating like a snake
Ring:
N monomers, N FENE bonds (terminal monomer bonded to first)
Closed loop topology → minimum energy shape is a circle
In pores: circular conformation is incompatible with pore-throat geometry
Must deform against bending energy κ to squeeze through → trapping!
For a semiflexible ring in a pore of diameter d:
Trapping occurs when d < 2R_g (ring diameter)
Escape requires thermal/active fluctuations to overcome bending energy barrier ~ κ
This explains the activity-induced shrinking of semiflexible rings (Fig. 7, right panel): higher Pe → more collisions with obstacles → stronger inward transverse fluctuations → ring compresses.
4.5 Radius of Gyration
```
R_g² = (1/N) Σᵢ |rᵢ − r_COM|²
     = (1/2N²) Σᵢ Σⱼ |rᵢ − rⱼ|²
```
For a 2D freely jointed chain in good solvent: R_g ~ N^ν, ν = 0.75 (Flory exponent in 2D).  
For a rod: R_g ~ N (ν = 1).  
For a Gaussian chain: R_g ~ N^{1/2}.
Activity-induced swelling: active fluctuations effectively increase the temperature T → larger R_g.  
Confinement-induced shrinking (ring): pore walls impose inward force → smaller R_g.
4.6 Gyration Tensor and Asphericity
The gyration tensor S is:
```
S_αβ = (1/N) Σᵢ (rᵢ,α − r_COM,α)(rᵢ,β − r_COM,β)
```
Eigenvalues: λ₁ ≤ λ₂ (in 2D)
```
R_g² = λ₁ + λ₂
Asphericity A = (λ₂ − λ₁)² / (λ₁ + λ₂)²
```
A = 0: isotropic (ring, globule)  
A = 1: fully extended (rod)
For a semiflexible ring, A should remain close to 0 (circular shape). For a semiflexible linear chain, A → 1. This is confirmed in Fig. 8.
---
5. Run Classification Algorithm
The displacement-ratio algorithm (Al Harraq et al. 2026, Note S5) classifies each frame:
```
δ(t) = Δr_net(t, t+w) / Δr_path(t, t+w)
```
where Δr_net = |r(t+w) − r(t)| and Δr_path = Σ|r(tᵢ₊₁) − r(tᵢ)| over window w.
δ is the inverse tortuosity: δ = 1 for perfectly straight motion, δ → 0 for stationary/circular motion.
Thresholds (validated against manually labeled E. coli trajectories):
δ < 0.6 → trapped
speed < v_thresh → trapped
otherwise → running
The trapping probability is then:
```
P_trapped = (number of trapped frames) / (total frames)
```
---
6. Directed Transport and Symmetry Breaking (Bhattacharyay 2011)
6.1 The Dimer Model
Bhattacharyay (2011) considers a dimer (two particles bound by a spring) with asymmetric damping:
```
ẍ₁ = −(1−β) ẋ₁ − α(x₁−x₂) + √(2(1−β)kBT) η₁ − ∂Φ/∂x₁
ẍ₂ = −ẋ₂ + α(x₁−x₂) + √(2kBT) η₂ − ∂Φ/∂x₂
```
Damping constants: (1−β) for particle 1, 1 for particle 2. The asymmetry β ≠ 0 breaks spatial symmetry.
6.2 Exact Result for COM Velocity
In the overdamped limit with hard-core collisions (fast collision time ≪ slow Langevin time):
```
V_COM = −β/(1−β) × √(kBTα/2π)
```
This is remarkable: the center of mass drifts at constant velocity in equilibrium — no external drive. The key physics:
Symmetry breaking: asymmetric damping creates unequal effective spring constants
The harmonic interaction is "effectively asymmetric" due to the damping
Hard-core collisions preserve momentum → no contribution to COM velocity
Second law is NOT violated: no energy extraction; entropy is constant (uniform motion)
6.3 Connection to Active Matter
This model shows that even without activity (Pe = 0), a symmetry-broken system in thermal equilibrium can undergo directed transport. In active matter terms, this is analogous to a particle with Pe > 0 in a symmetry-broken geometry — e.g., a chiral swimmer, or an ABP in a ratchet potential.
An open research question (relevant to your upcoming IISER Pune project on defect dynamics): Can activity replace the asymmetric damping? Specifically, does an active dimer with Pe > 0 and asymmetric swimming speeds on each particle show directed transport that scales differently with α and T than Bhattacharyay's result?
Predicted scaling for an active dimer:
```
V_COM ~ (v₀ β) / (1 + Pe²/Pe*²)^{1/2}    (dimensional estimate)
```
where Pe* = v₀/√(kBTα) is the crossover Péclet number. At Pe ≫ Pe*, activity dominates and V_COM saturates.
---
7. Phase Diagram — Physical Interpretation
The phase diagram (Fig. 9) separates two regimes in (L, Pe) space:
Trapped phase (small L, large Pe):
Particle pushes hard into walls with force ~ Pe × kBT/σ
Wall collision time τ_wall ~ L/v₀ ~ L/Pe → 0 as Pe↑
No time to reorient before hitting next wall
P_trapped → high
Run-and-tumble phase (large L, small Pe):
Free tumbling reorients particle before wall contact
L > ℓ_p = v₀/D_R = Pe/D_R → pore fits many run lengths
P_trapped → low
Phase boundary:
The transition roughly follows L* ~ ℓ_p = Pe (in reduced units with D_R = 1):
```
L*(Pe) ≈ Pe / D_R
```
At the boundary, the crossover length λ_c equals L* ≈ ℓ_p. This is the line where the Datta-group experimental results show the transition from loamy to silty behavior for E. coli.
---
8. Numerical Integration — Euler-Maruyama Scheme
The stochastic differential equation:
```
dx = f(x, θ) dt + √(2kBT) dW
```
is discretized as:
```
x(t+dt) = x(t) + f(x(t), θ(t)) × dt + √(2kBT dt) × N(0,1)
θ(t+dt) = θ(t) + √(2DR dt) × N(0,1)
```
Stability condition: dt < min(1/Pe, 1/DR, σ²/(2kBT))
In reduced units (kBT=1, σ=1, DR=1): dt < 0.5. We use dt = 10⁻³ (well within stability).
Recording interval: We save every 20 steps → dt_rec = 0.02 τ. This is much smaller than τ_R = 1/DR = 1 τ, so the MSD ballistic-to-diffusive crossover is well resolved.
---
9. Key Physical Parameters — Summary Table
Symbol	Meaning	Typical range	Reduced units
Pe	Péclet number	1–100	= v₀ (kBT=γ=σ=1)
D_R	Rotational diffusion	0.01–1 s⁻¹	= 1 (τ_R = 1)
τ_R	Persistence time	1–100 τ	1/D_R
ℓ_p	Persistence length	5–50 σ	Pe/D_R
L	Characteristic pore length	5–25 σ	measured
κ	Bending rigidity	0–1000 kBT	dimensionless
k_f	FENE spring constant	30–1000 kBT/σ²	dimensionless
P_trapped	Trapping probability	0–1	dimensionless
D_eff	Effective diffusivity	10⁻² – 1 σ²/τ	σ²/τ
λ_c	Crossover length	1–20 σ	σ
A	Asphericity	0–1	dimensionless
R_g	Radius of gyration	3–20 σ	σ
---
10. Suggestions for Extension (B.S.-M.S. Project Level)
10.1 Chemotaxis in Porous Media
Add a chemical field c(x,t) satisfying diffusion equation ∂c/∂t = D_c ∇²c − κ_deg c + S(x) (source at root). Particle tumbling rate λ(c) = λ₀ × f(c) where f is the chemotactic response function. This directly extends the Datta group's PRX Life 2026 results to quantify how the trade-off between confinement and chemical gradient determines recruitment efficiency.
10.2 Active Dimer with Asymmetric Activity
Implement Bhattacharyay's dimer with active forces: each particle has self-propulsion v₀_1 ≠ v₀_2. Measure V_COM as a function of (v₀_1 − v₀_2) and compare to the thermal result V ~ β√(αkBT). This is a new result not in the literature.
10.3 Defect Dynamics in Confined Active Nematics
Replace single ABPs with active nematics (head-tail symmetric pushers/pullers). Study ±1/2 defect pairs in porous confinement. This connects directly to your planned IISER Pune project on defect dynamics in active media under Prof. Hikkadi.
10.4 Predator-Prey in Porous Media
Implement two species: prey (ABP with Pe_prey) and predator (ABP with Pe_pred > Pe_prey + chemotaxis toward prey). Study how pore geometry modulates predator success rate. This connects to your Dr. Rumi De project on cooperative interactions.
