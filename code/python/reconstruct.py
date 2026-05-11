"""
Reconstruct marginal quantile functions and the trivariate copula density
from the published `*_coefficients_{normal,average}.csv` files.

These CSVs contain the polynomial coefficients of the distributional model
(see `Reconstruction.jl:108` in this repo). Columns x1..x1694 are the mutable
copula coefficients kappa, x1695..x1730 are the marginal quantile coefficients
xi (12 orders per measure, in alphabetical order: consum, income, wealth).

The basis is the orthonormal Legendre family on [0,1]:

    Q_o(u) = sqrt(2*o+1) * P_o(2*u - 1)

where P_o is the standard Legendre polynomial on [-1,1].

Reconstructions:

    quantile_at(date, measure, u)
        = Xi^{-1}_{m,t}(u) = sum_o xi[o, m, t] * Q_o(u)

    copula_density_at(date, u_c, u_y, u_w)
        = dC_t(u_c, u_y, u_w)
        = sum_{o_c, o_y, o_w} kappa[o_c, o_y, o_w, t] * Q_{o_c}(u_c) Q_{o_y}(u_y) Q_{o_w}(u_w)

CLI usage:

    python reconstruct.py PSID_coefficients_normal.csv quantile --date 2008-Q3 --measure consum --u 0.1,0.5,0.9
    python reconstruct.py PSID_coefficients_normal.csv copula   --date 2008-Q3 --u 0.5,0.5,0.5

Importable usage:

    from reconstruct import Reconstruction
    r = Reconstruction("PSID_coefficients_normal.csv")
    q = r.quantile_at("2008-Q3", "consum", [0.1, 0.5, 0.9])     # marginal quantile values
    d = r.copula_density_at("2008-Q3", 0.5, 0.5, 0.5)           # scalar density at middle
    grid = r.copula_density_grid("2008-Q3", n=30)               # 30x30x30 density values
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Model constants — match the Julia code (Reconstruction.jl / ModelPrep.jl)
# -----------------------------------------------------------------------------
GRID_COP = 12            # cop_part = 12^3 = 1728
GRID_PCF = 12            # 12 polynomial coefs per marginal
D = 3                    # consum, income, wealth
N_MUTABLE_COP = 1694     # 12^3 - 34
N_QUANTILE = D * GRID_PCF  # 36
MEASURES = ("consum", "income", "wealth")  # alphabetical, matches sort!(measures)


# -----------------------------------------------------------------------------
# Orthonormal Legendre basis on [0,1]
# -----------------------------------------------------------------------------
def legendre(m: int, x: np.ndarray) -> np.ndarray:
    """Standard Legendre polynomial P_m on [-1, 1]. Vectorized in x."""
    if m == 0:
        return np.ones_like(x, dtype=float)
    if m == 1:
        return np.asarray(x, dtype=float)
    p_prev_prev = np.ones_like(x, dtype=float)
    p_prev = np.asarray(x, dtype=float)
    for n in range(2, m + 1):
        p_curr = ((2 * n - 1) * x * p_prev - (n - 1) * p_prev_prev) / n
        p_prev_prev, p_prev = p_prev, p_curr
    return p_curr


def Q(m: int, u: np.ndarray) -> np.ndarray:
    """Orthonormal Legendre basis Q_m on [0, 1]."""
    u = np.asarray(u, dtype=float)
    return np.sqrt(2 * m + 1) * legendre(m, 2 * u - 1)


# -----------------------------------------------------------------------------
# Index mapping: (x1..x1694) -> mutable CartesianIndex in the 12^3 copula tensor
# -----------------------------------------------------------------------------
def _build_mutable_index_map() -> list[tuple[int, int, int]]:
    """
    Returns the ordered list of mutable 1-indexed CartesianIndices, matching
    Julia's `filter(!f, CartesianIndices((12,12,12)))` iteration.

    An entry (i, j, k) is "immutable" iff at least D-1 = 2 of its components == 1.
    Iteration is column-major (i varies fastest, then j, then k).
    """
    out = []
    for k in range(1, GRID_COP + 1):
        for j in range(1, GRID_COP + 1):
            for i in range(1, GRID_COP + 1):
                n_ones = (i == 1) + (j == 1) + (k == 1)
                if n_ones < D - 1:  # mutable
                    out.append((i, j, k))
    assert len(out) == N_MUTABLE_COP, f"expected {N_MUTABLE_COP}, got {len(out)}"
    return out


_MUTABLE_INDICES = _build_mutable_index_map()


def _unflatten_to_kappa(row: np.ndarray) -> np.ndarray:
    """
    Take the first 1694 entries of a row (mutable copula coefs) and place them
    into a (12, 12, 12) tensor at their (o_c, o_y, o_w) positions.

    The leading immutable entry kappa[0,0,0] is set to 1 (uniform marginal).
    All other immutable entries (rays where at least D-1 indices are 0 in
    polynomial-order space) are set to 0.
    """
    kappa = np.zeros((GRID_COP, GRID_COP, GRID_COP), dtype=float)
    for x_k, (i, j, k) in zip(row[:N_MUTABLE_COP], _MUTABLE_INDICES):
        kappa[i - 1, j - 1, k - 1] = x_k  # polynomial orders are 0-indexed
    kappa[0, 0, 0] = 1.0
    return kappa


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------
@dataclass
class Reconstruction:
    """Lazy reader for a *_coefficients_*.csv file."""

    csv_path: str | Path

    def __post_init__(self) -> None:
        self.df = pd.read_csv(self.csv_path)
        if "time" not in self.df.columns:
            raise ValueError("Expected a 'time' column.")
        n_data = self.df.shape[1] - 1
        expected = N_MUTABLE_COP + N_QUANTILE
        if n_data != expected:
            raise ValueError(
                f"{self.csv_path} has {n_data} data columns; expected {expected}. "
                "Is this a *_coefficients_*.csv file?"
            )
        self._dates = self.df["time"].astype(str).tolist()

    # --- helpers ---------------------------------------------------------------

    def available_dates(self) -> list[str]:
        return list(self._dates)

    def _row_for_date(self, date: str) -> np.ndarray:
        if date not in self._dates:
            raise KeyError(
                f"{date!r} not in file. First/last available: "
                f"{self._dates[0]!r} .. {self._dates[-1]!r}"
            )
        idx = self._dates.index(date)
        return self.df.iloc[idx, 1:].to_numpy(dtype=float)

    def _xi(self, date: str, measure: str) -> np.ndarray:
        """Return the 12 quantile polynomial coefficients for the given measure."""
        if measure not in MEASURES:
            raise ValueError(f"measure must be one of {MEASURES}, got {measure!r}")
        m_idx = MEASURES.index(measure)
        row = self._row_for_date(date)
        start = N_MUTABLE_COP + m_idx * GRID_PCF
        return row[start : start + GRID_PCF]

    def _kappa(self, date: str) -> np.ndarray:
        """Return the (12, 12, 12) copula coefficient tensor at the given date."""
        return _unflatten_to_kappa(self._row_for_date(date))

    # --- public methods -------------------------------------------------------

    def quantile_at(self, date: str, measure: str, u) -> np.ndarray:
        """
        Marginal quantile function Xi^{-1}_{m,t}(u) for one date and one measure,
        evaluated at quantile point(s) u in [0,1]. u may be a scalar or array.
        """
        xi = self._xi(date, measure)
        u = np.atleast_1d(np.asarray(u, dtype=float))
        # Build basis matrix [n_u, n_orders]
        basis = np.stack([Q(o, u) for o in range(GRID_PCF)], axis=-1)
        out = basis @ xi
        return out if out.ndim > 0 else float(out)

    def copula_density_at(self, date: str, u_c, u_y, u_w) -> float | np.ndarray:
        """
        Copula density dC_t(u_c, u_y, u_w). Scalar inputs return a scalar;
        equal-length 1D inputs return a 1D array (element-wise evaluation).
        """
        kappa = self._kappa(date)
        u_c, u_y, u_w = (np.atleast_1d(np.asarray(x, dtype=float)) for x in (u_c, u_y, u_w))
        if not (u_c.shape == u_y.shape == u_w.shape):
            raise ValueError("u_c, u_y, u_w must have the same shape")
        # basis[o, k] = Q_o(u[k]) — shape (12, n_points)
        B_c = np.stack([Q(o, u_c) for o in range(GRID_COP)], axis=0)
        B_y = np.stack([Q(o, u_y) for o in range(GRID_COP)], axis=0)
        B_w = np.stack([Q(o, u_w) for o in range(GRID_COP)], axis=0)
        # einsum over polynomial-order axes: kappa[a,b,c] * B_c[a,k] * B_y[b,k] * B_w[c,k] -> [k]
        out = np.einsum("abc,ak,bk,ck->k", kappa, B_c, B_y, B_w)
        return float(out[0]) if out.size == 1 else out

    def copula_density_grid(self, date: str, n: int = 30) -> np.ndarray:
        """
        Evaluate the copula density on a regular n x n x n grid of (u_c, u_y, u_w)
        in (0, 1). Useful for plotting heatmaps / iso-surfaces.

        Returns an ndarray of shape (n, n, n).
        """
        kappa = self._kappa(date)
        u = np.linspace(1e-6, 1 - 1e-6, n)
        B = np.stack([Q(o, u) for o in range(GRID_COP)], axis=0)  # (12, n)
        # kappa[a,b,c] * B[a,i] * B[b,j] * B[c,k] -> grid[i, j, k]
        return np.einsum("abc,ai,bj,ck->ijk", kappa, B, B, B)

    def marginal_grid(self, date: str, measure: str, n: int = 200) -> tuple[np.ndarray, np.ndarray]:
        """
        Evaluate the marginal quantile function on a regular grid u in (0, 1).
        Returns (u_grid, quantile_values).
        """
        u = np.linspace(1e-6, 1 - 1e-6, n)
        return u, self.quantile_at(date, measure, u)


# -----------------------------------------------------------------------------
# Sanity checks (run with `python -m pytest reconstruct.py` or as `__main__`)
# -----------------------------------------------------------------------------
def _self_check() -> None:
    """Basic invariants — does not require a CSV."""
    # 1) Legendre orthonormality on [0,1]: <Q_i, Q_j> = delta_ij (Monte Carlo, loose).
    rng = np.random.default_rng(0)
    u = rng.uniform(size=400_000)
    for i in range(4):
        for j in range(4):
            inner = float(np.mean(Q(i, u) * Q(j, u)))
            expected = 1.0 if i == j else 0.0
            assert abs(inner - expected) < 0.02, (i, j, inner)

    # 2) Index map: 1694 entries, none with >=2 ones, all unique.
    assert len(_MUTABLE_INDICES) == N_MUTABLE_COP
    assert len(set(_MUTABLE_INDICES)) == N_MUTABLE_COP
    for (i, j, k) in _MUTABLE_INDICES:
        assert ((i == 1) + (j == 1) + (k == 1)) < D - 1

    print("self-check OK")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def _parse_u_list(s: str) -> np.ndarray:
    return np.array([float(x) for x in s.split(",")])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("csv", help="Path to a *_coefficients_normal.csv or *_average.csv")
    p.add_argument(
        "mode",
        choices=("quantile", "copula", "grid", "selfcheck"),
        help="What to compute",
    )
    p.add_argument("--date", help="Quarter, e.g. 2008-Q3")
    p.add_argument("--measure", choices=MEASURES, help="Marginal (for `quantile`)")
    p.add_argument(
        "--u",
        type=_parse_u_list,
        default=np.array([0.1, 0.25, 0.5, 0.75, 0.9]),
        help="Comma-separated quantile point(s). For `copula`: must be 3 numbers",
    )
    p.add_argument("--n", type=int, default=30, help="Grid size for `grid` mode")
    args = p.parse_args()

    if args.mode == "selfcheck":
        _self_check()
        return

    r = Reconstruction(args.csv)

    if args.mode == "quantile":
        if not args.date or not args.measure:
            p.error("--date and --measure are required for `quantile`")
        q = r.quantile_at(args.date, args.measure, args.u)
        print("u\tquantile")
        for u_i, q_i in zip(args.u, np.atleast_1d(q)):
            print(f"{u_i:.4f}\t{q_i:.6g}")

    elif args.mode == "copula":
        if not args.date:
            p.error("--date is required for `copula`")
        if len(args.u) != 3:
            p.error("`copula` needs exactly 3 values in --u, e.g. 0.5,0.5,0.5")
        d = r.copula_density_at(args.date, args.u[0], args.u[1], args.u[2])
        print(f"dC({args.u[0]}, {args.u[1]}, {args.u[2]}) = {float(d):.6g}")

    elif args.mode == "grid":
        if not args.date:
            p.error("--date is required for `grid`")
        g = r.copula_density_grid(args.date, n=args.n)
        out = Path(args.csv).with_suffix("").name + f"_{args.date}_grid{args.n}.npy"
        np.save(out, g)
        print(f"saved {out}  shape={g.shape}  min={g.min():.4g}  max={g.max():.4g}")


if __name__ == "__main__":
    main()
