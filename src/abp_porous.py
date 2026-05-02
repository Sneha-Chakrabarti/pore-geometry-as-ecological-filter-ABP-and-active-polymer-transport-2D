
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:45:23 2026

@author: chakr
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
from pathlib import Path
import matplotlib.pyplot as plt
@dataclass
class PorousMedia:
    """Randomly placed, non-overlapping circular obstacles in a 2D periodic box."""
    box_size: float
    n_obstacles: int
    r_obs: float
    seed: int = 42
    centers: np.ndarray = field(init=False)

    def __post_init__(self):
        rng = np.random.default_rng(self.seed)
        centers = []
        attempts = 0
        max_attempts = 100_000
        while len(centers) < self.n_obstacles and attempts < max_attempts:
            c = rng.uniform(0, self.box_size, 2)
            if self._no_overlap(c, centers):
                centers.append(c)
            attempts += 1
        if len(centers) < self.n_obstacles:
            raise ValueError(
                f"Could only place {len(centers)}/{self.n_obstacles} obstacles. "
                "Reduce obstacle density or radius."
            )
        self.centers = np.array(centers)

    def _no_overlap(self, c, existing):
        if not existing:
            return True
        diffs = np.array(existing) - c
        dists = np.linalg.norm(diffs, axis=1)
        return np.all(dists > 2 * self.r_obs * 1.05)   # 5 % clearance

    def characteristic_pore_length(self, n_chords: int = 5000, rng_seed: int = 0):
        """
        Estimate characteristic pore length L by random-chord sampling.
        Shoot random line segments; record intersection-free lengths.
        Fits exponential: P(l) = (1/L) exp(-l/L)  (Torquato & Lu, Phys Rev E 1993).
        """
        rng = np.random.default_rng(rng_seed)
        chord_lengths = []
        step = 0.5 * self.r_obs
        for _ in range(n_chords):
            origin = rng.uniform(0, self.box_size, 2)
            angle  = rng.uniform(0, 2 * np.pi)
            direction = np.array([np.cos(angle), np.sin(angle)])
            # walk until hitting an obstacle or crossing the box
            t = 0.0
            hit = False
            while t < self.box_size * np.sqrt(2):
                pos = origin + t * direction
                # periodic wrap
                pos_w = pos % self.box_size
                dists = np.linalg.norm(self.centers - pos_w, axis=1)
                if np.any(dists < self.r_obs):
                    hit = True
                    break
                t += step
            if hit and t > step:
                chord_lengths.append(t)
        if not chord_lengths:
            return np.nan
        chord_lengths = np.array(chord_lengths)
        # MLE estimate for exponential distribution
        L = np.mean(chord_lengths)
        return L


@dataclass
class ABPSimulation:
    """
    Single active Brownian particle in 2D porous media.

    Parameters
    ----------
    Pe          : Péclet number  Pe = Fa * sigma / (kB*T)
    D_R         : rotational diffusion coefficient  [1/tau]
    dt          : integration timestep              [tau]
    media       : PorousMedia instance
    sigma       : particle diameter                 [sigma]  (=1 by default)
    """
    Pe: float
    D_R: float
    media: PorousMedia
    dt: float = 1e-3
    sigma: float = 1.0
    seed: int = 0

    # derived
    v0: float = field(init=False)
    tau_R: float = field(init=False)

    def __post_init__(self):
        self.v0  = self.Pe          # in reduced units where kBT=1, gamma=1, sigma=1
        self.tau_R = 1.0 / self.D_R

    def _find_free_position(self, rng):
        """Place particle away from obstacles."""
        for _ in range(100_000):
            pos = rng.uniform(0, self.media.box_size, 2)
            dists = np.linalg.norm(self.media.centers - pos, axis=1)
            if np.all(dists > (self.media.r_obs + 0.5 * self.sigma)):
                return pos
        raise RuntimeError("Cannot find free starting position.")

    def run(self, n_steps: int, record_every: int = 10):
        """
        Integrate overdamped Langevin + rotational diffusion.

        Returns
        -------
        positions  : (N_rec, 2)  — unwrapped positions for MSD
        wrapped    : (N_rec, 2)  — wrapped positions (inside box)
        thetas     : (N_rec,)    — orientation angle
        """
        rng = np.random.default_rng(self.seed)
        pos   = self._find_free_position(rng)
        theta = rng.uniform(0, 2 * np.pi)

        positions = []
        wrapped   = []
        thetas    = []
        # keep running sum for unwrapped trajectory
        cumulative = pos.copy()

        sqrt2DR_dt = np.sqrt(2 * self.D_R * self.dt)
        sqrt2kBT_dt = np.sqrt(2.0 * self.dt)   # kBT=1, gamma=1

        prev_wrapped = pos.copy()

        for step in range(n_steps):
            # ── Active force direction ──
            e = np.array([np.cos(theta), np.sin(theta)])

            # ── Thermal translational noise ──
            xi = rng.standard_normal(2)

            # ── Obstacle repulsion (WCA-like: steep soft wall) ──
            f_obs = np.zeros(2)
            diffs = pos - self.media.centers      # not periodic for simplicity (box >> obs)
            dists = np.linalg.norm(diffs, axis=1)
            contact = self.media.r_obs + 0.5 * self.sigma
            for i, d in enumerate(dists):
                if 0 < d < contact * 1.12:        # within WCA cutoff
                    overlap = contact / d
                    f_obs += 4.0 * (12 * overlap**12 - 6 * overlap**6) / d**2 * diffs[i]

            # ── Euler-Maruyama step ──
            dp = (self.v0 * e + f_obs) * self.dt + sqrt2kBT_dt * xi

            new_pos = (pos + dp) % self.media.box_size

            # track unwrapped displacement
            delta = new_pos - prev_wrapped
            # correct for PBC jump
            delta = delta - self.media.box_size * np.round(delta / self.media.box_size)
            cumulative += delta
            prev_wrapped = new_pos.copy()
            pos = new_pos

            # ── Rotational diffusion ──
            theta += sqrt2DR_dt * rng.standard_normal()

            if step % record_every == 0:
                positions.append(cumulative.copy())
                wrapped.append(pos.copy())
                thetas.append(theta)

        return (
            np.array(positions),
            np.array(wrapped),
            np.array(thetas)
        )


# ──────────────────────────────────────────────────────────────
# Analysis utilities
# ──────────────────────────────────────────────────────────────

def compute_msd(positions: np.ndarray, max_lag_fraction: float = 0.3):
    """
    Compute ensemble/time-averaged MSD from a single trajectory.

    MSD(tau) = < |r(t+tau) - r(t)|^2 >_t

    Returns lag indices and MSD values.
    """
    N = len(positions)
    max_lag = max(1, int(N * max_lag_fraction))
    lags = np.arange(1, max_lag + 1)
    msd  = np.zeros(len(lags))
    for k, lag in enumerate(lags):
        diff = positions[lag:] - positions[:-lag]
        msd[k] = np.mean(np.sum(diff**2, axis=1))
    return lags, msd


def classify_frames(positions: np.ndarray, window: int = 20, threshold: float = 0.6,
                    v_thresh: float = 0.05, dt_rec: float = 0.01):
    """
    Classify each frame as 'running' or 'trapped'.

    Mirrors the algorithm of Al Harraq et al. PRX Life 2026:
      displacement_ratio δ(t) = net_disp / path_length  over window w
      trapped if δ < threshold  OR  instantaneous speed < v_thresh
    """
    N = len(positions)
    labels = np.zeros(N, dtype=int)   # 0 = running, 1 = trapped

    for i in range(N):
        i_start = max(0, i - window // 2)
        i_end   = min(N - 1, i + window // 2)
        segment = positions[i_start:i_end + 1]
        if len(segment) < 2:
            continue
        net_disp   = np.linalg.norm(segment[-1] - segment[0])
        path_len   = np.sum(np.linalg.norm(np.diff(segment, axis=0), axis=1))
        if path_len < 1e-10:
            labels[i] = 1
            continue
        delta = net_disp / path_len
        # instantaneous speed at frame i
        if i > 0:
            speed = np.linalg.norm(positions[i] - positions[i-1]) / dt_rec
        else:
            speed = v_thresh + 1
        if delta < threshold or speed < v_thresh:
            labels[i] = 1
    return labels


def trapping_probability(labels: np.ndarray) -> float:
    """Fraction of frames spent trapped."""
    return np.mean(labels)


def effective_diffusivity(lags: np.ndarray, msd: np.ndarray,
                          dt_rec: float, fraction: float = 0.3) -> float:
    """
    Extract D_eff from the diffusive (linear) tail of the MSD.
    D_eff = MSD / (4 * tau)   [2D]
    Uses the last `fraction` of lags for the linear fit.
    """
    n = max(2, int(len(lags) * fraction))
    t = lags[-n:] * dt_rec
    m = msd[-n:]
    coeffs = np.polyfit(t, m, 1)
    return coeffs[0] / 4.0


def crossover_length(lags: np.ndarray, msd: np.ndarray, dt_rec: float):
    """
    Find the ballistic-to-diffusive crossover by segmented log-log fit.
    Returns crossover length lambda_c.
    """
    log_t = np.log10(lags * dt_rec + 1e-15)
    log_m = np.log10(msd + 1e-15)
    N = len(log_t)
    best_err = np.inf
    best_k   = N // 2
    for k in range(3, N - 3):
        # ballistic segment [0:k]
        b1 = np.polyfit(log_t[:k], log_m[:k], 1)
        # diffusive segment [k:]
        b2 = np.polyfit(log_t[k:], log_m[k:], 1)
        err = (np.sum((np.polyval(b1, log_t[:k]) - log_m[:k])**2) +
               np.sum((np.polyval(b2, log_t[k:]) - log_m[k:])**2))
        if err < best_err:
            best_err = err
            best_k   = k
    t_c = (lags[best_k] * dt_rec)
    m_c = msd[best_k]
    lambda_c = np.sqrt(m_c)
    return lambda_c, t_c
