
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:47:54 2026

@author: chakr
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Literal
from abp_porous import PorousMedia


@dataclass
class ActivePolymer:
    """
    Active Brownian polymer (linear or ring) in 2D porous media.

    Parameters (all in reduced units: sigma=kBT=gamma=1)
    ----------
    N        : number of monomers
    topology : 'linear' or 'ring'
    Pe       : Péclet number per monomer
    kf       : FENE spring constant
    kappa    : bending rigidity
    D_R      : rotational diffusion coefficient
    media    : PorousMedia instance
    dt       : timestep
    """
    N: int
    topology: Literal['linear', 'ring']
    Pe: float
    kf: float
    kappa: float
    D_R: float
    media: PorousMedia
    dt: float = 5e-4
    r_max: float = 1.5        # FENE max extension
    seed: int = 0

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed)
        self._init_chain()

    def _init_chain(self):
        """Place chain/ring in a pore region."""
        bs = self.media.box_size
        for _ in range(10_000):
            start = self.rng.uniform(bs * 0.1, bs * 0.9, 2)
            positions = [start]
            ok = True
            for _ in range(self.N - 1):
                angle = self.rng.uniform(0, 2*np.pi)
                nxt = positions[-1] + np.array([np.cos(angle), np.sin(angle)]) * 0.95
                nxt = nxt % bs
                dists = np.linalg.norm(self.media.centers - nxt, axis=1)
                if np.any(dists < self.media.r_obs + 0.5):
                    ok = False
                    break
                positions.append(nxt)
            if ok:
                self.pos = np.array(positions, dtype=float)
                self.theta = self.rng.uniform(0, 2*np.pi, self.N)
                return
        raise RuntimeError("Cannot initialise chain — reduce density.")

    # ── Bond list ───────────────────────────────────────────────
    def _bond_pairs(self):
        pairs = [(i, i+1) for i in range(self.N-1)]
        if self.topology == 'ring':
            pairs.append((self.N-1, 0))
        return pairs

    # ── Forces ──────────────────────────────────────────────────

    def _fene_force(self):
        """FENE spring between bonded monomers."""
        f = np.zeros((self.N, 2))
        bs = self.media.box_size
        for i, j in self._bond_pairs():
            dr = self.pos[j] - self.pos[i]
            dr -= bs * np.round(dr / bs)   # MIC
            r  = np.linalg.norm(dr)
            r  = max(r, 1e-8)
            if r >= self.r_max:
                r = self.r_max * 0.999
            ratio = r / self.r_max
            k_eff = self.kf / (1 - ratio**2)
            fij   = -k_eff * dr
            f[i] -= fij
            f[j] += fij
        return f

    def _bending_force(self):
        """Kratky-Porod bending rigidity."""
        f   = np.zeros((self.N, 2))
        if self.kappa == 0.0:
            return f
        bs  = self.media.box_size
        pairs = self._bond_pairs()
        n_bonds = len(pairs)
        for k in range(n_bonds - 1):
            i, j = pairs[k]
            _, l  = pairs[k+1]
            b1 = self.pos[j] - self.pos[i]; b1 -= bs*np.round(b1/bs)
            b2 = self.pos[l] - self.pos[j]; b2 -= bs*np.round(b2/bs)
            r1 = np.linalg.norm(b1); r2 = np.linalg.norm(b2)
            if r1 < 1e-8 or r2 < 1e-8:
                continue
            b1h = b1/r1; b2h = b2/r2
            cos_phi = np.clip(np.dot(b1h, b2h), -1, 1)
            # gradient of kappa*(1-cos_phi) w.r.t. positions
            fi = -self.kappa / r1 * (b2h - cos_phi * b1h)
            fl =  self.kappa / r2 * (b1h - cos_phi * b2h)
            fj = -(fi + fl)
            f[i] += fi; f[j] += fj; f[l] += fl
        return f

    def _wca_monomer_force(self):
        """WCA repulsion between non-bonded monomers."""
        f   = np.zeros((self.N, 2))
        bs  = self.media.box_size
        cutoff = 2**(1/6)
        for i in range(self.N):
            for j in range(i+2, self.N):
                if self.topology == 'ring' and i == 0 and j == self.N-1:
                    continue   # bonded pair
                dr = self.pos[j] - self.pos[i]; dr -= bs*np.round(dr/bs)
                r  = np.linalg.norm(dr)
                if r < cutoff and r > 1e-8:
                    sr6  = (1.0/r)**6
                    fmag = 4*(12*sr6**2/r - 6*sr6/r) / r
                    fij  = fmag * dr
                    f[i] -= fij; f[j] += fij
        return f

    def _obstacle_force(self):
        """WCA repulsion from obstacle surfaces."""
        f   = np.zeros((self.N, 2))
        bs  = self.media.box_size
        for i in range(self.N):
            for oc in self.media.centers:
                dr   = self.pos[i] - oc
                dr  -= bs * np.round(dr/bs)
                r    = np.linalg.norm(dr)
                sig  = self.media.r_obs + 0.5
                cutoff = sig * 2**(1/6)
                if r < cutoff and r > 1e-8:
                    sr6  = (sig/r)**6
                    fmag = 4*(12*sr6**2/r - 6*sr6/r) / r
                    f[i] += fmag * dr
        return f

    def _active_force(self):
        """Active propulsion along orientation vector."""
        e = np.stack([np.cos(self.theta), np.sin(self.theta)], axis=1)
        return self.Pe * e

    # ── Integration ─────────────────────────────────────────────

    def step(self):
        bs = self.media.box_size
        sqrt2dt = np.sqrt(2 * self.dt)
        sqrt2DRdt = np.sqrt(2 * self.D_R * self.dt)

        total_f = (self._fene_force() +
                   self._bending_force() +
                   self._wca_monomer_force() +
                   self._obstacle_force() +
                   self._active_force())

        xi  = self.rng.standard_normal((self.N, 2))
        dpos = total_f * self.dt + sqrt2dt * xi
        self.pos = (self.pos + dpos) % bs
        self.theta += sqrt2DRdt * self.rng.standard_normal(self.N)

    def run(self, n_steps: int, record_every: int = 100):
        """Run simulation; record COM position (unwrapped) and Rg."""
        bs = self.media.box_size
        com_list, rg_list = [], []
        cum_com = self._com().copy()
        prev_com = self._com().copy()

        for s in range(n_steps):
            self.step()
            if s % record_every == 0:
                com_w = self._com()
                delta = com_w - prev_com
                delta -= bs * np.round(delta / bs)
                cum_com += delta
                prev_com = com_w.copy()
                com_list.append(cum_com.copy())
                rg_list.append(self._rg())
        return np.array(com_list), np.array(rg_list)

    # ── Observables ─────────────────────────────────────────────

    def _com(self):
        bs = self.media.box_size
        # unwrap relative to first monomer using MIC
        ref = self.pos[0]
        dp  = self.pos - ref
        dp -= bs * np.round(dp / bs)
        return (ref + np.mean(dp, axis=0)) % bs

    def _rg(self):
        """Radius of gyration."""
        bs = self.media.box_size
        ref = self.pos[0]
        dp  = self.pos - ref; dp -= bs*np.round(dp/bs)
        unwr = ref + dp
        com  = np.mean(unwr, axis=0)
        diff = unwr - com
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    def gyration_tensor(self):
        """2x2 gyration tensor — eigenvalues give asphericity."""
        bs   = self.media.box_size
        ref  = self.pos[0]
        dp   = self.pos - ref; dp -= bs*np.round(dp/bs)
        unwr = ref + dp
        com  = np.mean(unwr, axis=0)
        diff = unwr - com
        S    = diff.T @ diff / self.N
        return S

    def asphericity(self):
        """
        A = (lambda2 - lambda1)^2 / (lambda1+lambda2)^2
        A=0: circular; A=1: fully extended rod.
        """
        S  = self.gyration_tensor()
        ev = np.linalg.eigvalsh(S)
        l1, l2 = sorted(ev)
        denom = (l1 + l2)**2
        if denom < 1e-12:
            return 0.0
        return (l2 - l1)**2 / denom


# ── Batch runners ────────────────────────────────────────────────────────────

def run_msd_sweep(Pe_values, topology, kf, kappa, D_R, media,
                  n_steps=200_000, record_every=50, dt=5e-4, N=20, seed_base=0):
    """
    For each Pe, run one trajectory and return (lags, msd) dict.
    """
    from abp_porous import compute_msd
    results = {}
    dt_rec  = dt * record_every
    for Pe in Pe_values:
        poly = ActivePolymer(N=N, topology=topology, Pe=Pe, kf=kf,
                              kappa=kappa, D_R=D_R, media=media,
                              dt=dt, seed=seed_base)
        coms, _ = poly.run(n_steps, record_every)
        lags, msd = compute_msd(coms, max_lag_fraction=0.25)
        results[Pe] = (lags * dt_rec, msd)
    return results


def run_rg_sweep(Pe_values, topology, kf, kappa, D_R, media,
                 n_steps=200_000, record_every=50, dt=5e-4, N=20, seed_base=0):
    """Return mean Rg for each Pe."""
    results = {}
    for Pe in Pe_values:
        poly = ActivePolymer(N=N, topology=topology, Pe=Pe, kf=kf,
                              kappa=kappa, D_R=D_R, media=media,
                              dt=dt, seed=seed_base)
        _, rgs = poly.run(n_steps, record_every)
        results[Pe] = float(np.mean(rgs[len(rgs)//2:]))   # discard equilibration
    return results
