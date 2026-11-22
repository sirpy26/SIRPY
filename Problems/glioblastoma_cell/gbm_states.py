#!/usr/bin/env python3
"""Glioblastoma cell states under temozolomide: what the barcodes can and
cannot tell. A re-analysis, with SIRPY, of the public data of

  Larsson I, et al. Modeling glioblastoma heterogeneity as a dynamic network
  of cell states. Mol Syst Biol 17:e10105 (2021). Data and code:
  github.com/idalarsson/STAG (file expressionMatrices_20210622.mat).

    python3 gbm_states.py          # reads the two CSV files next to this
                                   # script and writes everything to out_gbm/

Three parts:
 1. Fates over 8 days (no differential equations): for one cell in each
    state at day 8, how many descendants in each state at day 16, under
    control, temozolomide (TMZ) and BMP4. Poisson likelihood with the
    authors' sampling fractions, extra variability (overdispersion) taken
    into account; exact profile-likelihood 95% intervals.
 2. Rates (SIRPY): the authors' model dX/dt = A X, fitted by maximum
    likelihood. For each of the 30 possible switches: does the data
    establish it (likelihood-ratio test against 'this switch is zero')?
 3. Treatment order (SIRPY): cells after 16 days for the two drugs in either
    order, built from the 8-day exposures that were measured, with
    approximate 95% intervals.

How SIRPY is used: every reported fit, test and prediction is computed from
SIRPY solutions of the model, and every solution must pass SIRPY's own check
('solved+verified'); derivatives come from the same solutions. Because this
model is linear, a closed-form matrix exponential exists; it is used only to
search for the best parameters quickly, and SIRPY then recomputes and
confirms every reported value.

Every number in the printout and on the figures is computed in this run.
"""
import json
import math
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

# the optimiser's trial steps can overflow harmlessly; keep the console clean
warnings.filterwarnings("ignore", category=RuntimeWarning)

import numpy as np
from scipy.linalg import expm, logm
from scipy.optimize import minimize

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

import sirpy  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out_gbm")
os.makedirs(OUT, exist_ok=True)
K = 6
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:7.1f} s]", *a, flush=True)


# ============================================================== data
def read_cells(path):
    import csv
    rows = list(csv.DictReader(open(path)))
    return [(r["barcode"], int(r["day"]), int(r["treatment"]), int(r["state"]))
            for r in rows]


def count_matrix(cells, treatment, day, barcodes):
    idx = {b: i for i, b in enumerate(barcodes)}
    X = np.zeros((K, len(barcodes)))
    for b, d, t, s in cells:
        if t == treatment and d == day:
            X[s - 1, idx[b]] += 1
    return X


def load_data():
    tc = read_cells(os.path.join(HERE, "U3065MG_timecourse_cells.csv"))
    tr = read_cells(os.path.join(HERE, "U3065MG_treatment_cells.csv"))
    out = {}
    # time course (untreated): samples at days 7, 14, 21 (coded 1, 2, 3);
    # authors' sampling fractions 0.8, 0.8, 1.0
    bcs = sorted({c[0] for c in tc})
    X = [count_matrix(tc, 1, d, bcs) for d in (1, 2, 3)]
    f = [0.8, 0.8, 1.0]
    out["untreated"] = [(X[a], X[b], (1 - f[a]) / f[a] * f[b], 7.0)
                        for a, b in ((0, 1), (1, 2))]
    # treatment experiment: day 8 (coded 0, shared) to day 16 (coded 8);
    # authors' fractions 0.7 and 1/3; treatment 1 control, 2 BMP4, 3 TMZ
    bcs = sorted({c[0] for c in tr})
    Xb = count_matrix(tr, 1, 0, bcs)
    c2 = (1 - 0.7) / 0.7 * (1 / 3)
    for t, name in ((1, "control"), (2, "BMP4"), (3, "TMZ")):
        out[name] = [(Xb, count_matrix(tr, t, 8, bcs), c2, 8.0)]
    # keep only lineages seen in the 'before' sample (others carry no
    # information about the dynamics in this model)
    for name, pairs in out.items():
        new = []
        for Xb_, Xa_, c2_, dt_ in pairs:
            keep = Xb_.sum(0) > 0
            new.append((Xb_[:, keep], Xa_[:, keep], c2_, dt_))
        out[name] = new
    x_day8 = Xb.sum(1) * (1 - 0.7) / 0.7       # kept day-8 population by state
    return out, x_day8


# ============================================================== part 1
class Fates:
    """8-day (or 7-day) fate matrices P >= 0, Poisson likelihood
    after ~ Poisson(c2 P before); one P per condition."""

    def __init__(self, conds):
        self.conds = conds                      # list of lists of pairs
        self.n = len(conds) * K * K

    def split(self, v):
        return [v[i * K * K:(i + 1) * K * K].reshape(K, K)
                for i in range(len(self.conds))]

    def f(self, v, a=None, lam=0.0):
        val, grad = 0.0, []
        for P, pairs in zip(self.split(v), self.conds):
            g = np.zeros((K, K))
            for xb, xa, c2, _ in pairs:
                mu = np.maximum(c2 * P @ xb, 1e-300)
                val += float(np.sum(mu - xa * np.log(mu)))
                g += c2 * (1.0 - xa / mu) @ xb.T
            grad.append(g.ravel())
        grad = np.concatenate(grad)
        if a is not None:
            val -= lam * float(a @ v)
            grad = grad - lam * a
        return val, grad

    def solve(self, a=None, lam=0.0, v0=None):
        v0 = np.full(self.n, 0.5) if v0 is None else np.maximum(v0, 1e-8)
        r = minimize(self.f, v0, args=(a, lam), jac=True, method="L-BFGS-B",
                     bounds=[(0.0, None)] * self.n,
                     options={"maxiter": 50000, "ftol": 1e-15,
                              "gtol": 1e-10, "maxcor": 30})
        return r.x

    def phi(self, v):
        s, n = 0.0, 0
        for P, pairs in zip(self.split(v), self.conds):
            for xb, xa, c2, _ in pairs:
                mu = c2 * P @ xb
                ok = mu > 1e-9
                s += float(np.sum((xa[ok] - mu[ok]) ** 2 / mu[ok]))
                n += int(ok.sum())
        return s / (n - self.n)

    def interval(self, a, v_hat, phi):
        """Exact profile-likelihood 95% interval for <a, P> (convex problem),
        traced through Lagrange multipliers."""
        nll_hat = self.f(v_hat)[0]
        target = 3.8415 * phi / 2.0

        def side(sign):
            lo_l, hi_l, v = 0.0, sign * 1.0, v_hat
            for _ in range(60):
                v = self.solve(a, hi_l, v)
                if self.f(v)[0] - nll_hat > target:
                    break
                lo_l, hi_l = hi_l, hi_l * 2
            else:
                return None
            for _ in range(50):
                m = 0.5 * (lo_l + hi_l)
                v = self.solve(a, m, v)
                if self.f(v)[0] - nll_hat > target:
                    hi_l = m
                else:
                    lo_l = m
            return float(a @ self.solve(a, lo_l, v))

        return float(a @ v_hat), side(-1.0), side(+1.0)


# ============================================================== SIRPY model
SYSTEM = {f"Dy{i+1}": " + ".join(f"a{i}_{j}*y{j+1}" for j in range(K))
          for i in range(K)}
SOLVES = {"solves": 0, "verified": 0}
NQ = 40                                           # Gauss-Legendre nodes


class Propagator:
    """exp(sA) for every s in [0, Tmax], from K SIRPY solves (one per
    starting state); each must come back 'solved+verified'."""

    def __init__(self, A, Tmax):
        self.A, self.Tmax = np.asarray(A, float), float(Tmax)
        params = {f"a{i}_{j}": float(self.A[i, j])
                  for i in range(K) for j in range(K)}
        self.sol = []
        for c in range(K):
            ic = {f"y{i+1}": (1.0 if i == c else 0.0) for i in range(K)}
            r = sirpy.solve_ivp(system=SYSTEM, x0=0, x1=self.Tmax, ic=ic,
                                params=params, verbose=False)
            SOLVES["solves"] += 1
            if r.get("status") != "solved+verified":
                raise RuntimeError("a SIRPY solve was not verified: "
                                   f"{r.get('status')}")
            SOLVES["verified"] += 1
            self.sol.append(r)
        self._c = {}

    def Phi(self, s):
        key = round(float(s), 14)
        if key not in self._c:
            M = np.empty((K, K))
            for c, r in enumerate(self.sol):
                for i in range(K):
                    M[i, c] = r["evaluate"](f"y{i+1}", float(s))
            self._c[key] = M
        return self._c[key]

    def nodes(self, t):
        x, w = np.polynomial.legendre.leggauss(NQ)
        return 0.5 * t * (x + 1.0), 0.5 * t * w

    def grad_A(self, GP, t):
        """d/dA <GP, exp(tA)> = int_0^t Phi(s)^T GP Phi(t-s)^T ds."""
        out = np.zeros((K, K))
        for s, w in zip(*self.nodes(t)):
            out += w * self.Phi(s).T @ GP @ self.Phi(t - s).T
        return out

    def dP_dA(self, t):
        """d exp(tA) / dA_ij for all i, j: array [i, j, :, :]."""
        D = np.zeros((K, K, K, K))
        for s, w in zip(*self.nodes(t)):
            D += w * np.einsum("ai,jb->ijab", self.Phi(s), self.Phi(t - s))
        return D


# parameters: g_j (net growth of state j) and s_ij (switch j -> i, >= 0)
OFF = [(i, j) for i in range(K) for j in range(K) if i != j]
NP = K + len(OFF)
LAB = [f"g{j+1}" for j in range(K)] + [f"s{j+1}->{i+1}" for i, j in OFF]
DA = np.zeros((NP, K, K))
for j in range(K):
    DA[j, j, j] = 1.0
for k, (i, j) in enumerate(OFF):
    DA[K + k, i, j] = 1.0
    DA[K + k, j, j] = -1.0


def to_A(th):
    A = np.einsum("p,pij->ij", np.asarray(th, float), DA)
    return A


def from_A(A):
    return np.array([A[:, j].sum() for j in range(K)]
                    + [A[i, j] for i, j in OFF])


class RateLik:
    def __init__(self, pairs):
        self.pairs = pairs
        self.Tmax = max(p[3] for p in pairs)

    def evaluate(self, th, fisher=True):
        pr = Propagator(to_A(th), self.Tmax)
        nll, gA = 0.0, np.zeros((K, K))
        F = np.zeros((NP, NP))
        pear, nobs = 0.0, 0
        for xb, xa, c2, dt in self.pairs:
            mu = c2 * pr.Phi(dt) @ xb
            if np.any(mu < -1e-9):
                return np.inf, None, None, None, pr
            mu = np.maximum(mu, 1e-300)
            nll += float(np.sum(mu - xa * np.log(mu)))
            gA += pr.grad_A(c2 * (1.0 - xa / mu) @ xb.T, dt)
            ok = mu > 1e-12
            pear += float(np.sum((xa[ok] - mu[ok]) ** 2 / mu[ok]))
            nobs += int(ok.sum())
            if fisher:
                dP = np.einsum("pij,ijab->pab", DA, pr.dP_dA(dt))
                J = (c2 * np.einsum("pab,bn->pan", dP, xb)).reshape(NP, -1)
                F += (J / mu.ravel()) @ J.T
        g = np.einsum("pij,ij->p", DA, gA)
        return nll, g, F, (pear, nobs), pr


def _qp(gf, H, lo, x):
    """min g.d + d.H.d/2 subject to x + d >= lo (H positive definite), solved
    exactly as a bounded least-squares problem: with H = L L^T this is
    min |L^T d + L^-1 g|^2 / 2 (bounded-variable least squares)."""
    from scipy.linalg import cholesky, solve_triangular
    from scipy.optimize import lsq_linear
    L = cholesky(H, lower=True)
    b = -solve_triangular(L, gf, lower=True)
    lb = np.where(np.isfinite(lo), lo - x, -np.inf)
    d = lsq_linear(L.T, b, bounds=(lb, np.full(len(x), np.inf)),
                   method="bvls", tol=1e-14, lsmr_tol=None).x
    return d, -float(gf @ d + 0.5 * d @ H @ d)


def lm_fit(lik, th0, fixed=None, gap=1e-6, maxit=300):
    """Bound-constrained Levenberg-Marquardt (Fisher scoring). Stops when the
    undamped model step predicts a decrease of the negative log-likelihood
    below `gap`; unlike a gradient test, this also covers flat directions."""
    fixed = dict(fixed or {})
    th = np.array(th0, float)
    for k, v in fixed.items():
        th[k] = v
    th[K:] = np.maximum(th[K:], 0.0)
    free = np.array([k for k in range(NP) if k not in fixed])
    lo = np.array([-np.inf if k < K else 0.0 for k in free])
    f, g, F, pe, _ = lik.evaluate(th)
    lam, pred0 = 1e-3, np.inf
    for _ in range(maxit):
        x, gf, Ff = th[free], g[free], F[np.ix_(free, free)]
        ridge = 1e-12 * np.max(np.diag(Ff)) * np.eye(len(free))
        _, pred0 = _qp(gf, Ff + ridge, lo, x)
        if pred0 < gap:
            break
        moved = False
        for _ in range(40):
            d, pred = _qp(gf, Ff + lam * np.diag(np.diag(Ff)) + ridge, lo, x)
            cand = th.copy()
            cand[free] = x + d
            cand[K:] = np.maximum(cand[K:], 0.0)
            fc, gc, Fc, pc, _ = lik.evaluate(cand)
            if np.isfinite(fc) and fc < f:
                rho = (f - fc) / max(pred, 1e-300)
                th, f, g, F, pe = cand, fc, gc, Fc, pc
                lam = min(max(lam * (0.2 if rho > 0.75 else 1.0 if rho > 0.25 else 3.0), 1e-12), 1e12)
                moved = True
                break
            lam = min(lam * 4.0, 1e12)
            if lam >= 1e12:
                break
        if not moved:
            break
    x = th[free]
    pg = g[free].copy()
    pg[(x <= lo + 1e-14) & (g[free] > 0)] = 0.0
    sd = np.sqrt(np.maximum(np.diag(F[np.ix_(free, free)]), 1e-12))
    return {"theta": th, "nll": f, "grad": g, "F": F, "pearson": pe,
            "score": float(np.max(np.abs(pg) / sd)), "gap": float(pred0)}


def start_from_fates(P, dt):
    try:
        L = logm(P)
        if np.max(np.abs(np.imag(L))) > 1e-8 * max(1.0, np.max(np.abs(L))) \
                or not np.all(np.isfinite(L)):
            raise ValueError
        A = np.real(L) / dt
    except Exception:
        A = (P - np.eye(K)) / dt
    off = ~np.eye(K, dtype=bool)
    A[off] = np.maximum(A[off], 1e-3)
    th = from_A(A)
    return th


def cf_fit(pairs, th0, fixed=None):
    """L-BFGS-B search with the closed-form matrix exponential, which exists
    for this linear model. Returns the parameter vector found."""
    from scipy.linalg import expm_frechet
    fixed = dict(fixed or {})
    base = np.array(th0, float)
    for k, v in fixed.items():
        base[k] = v
    free = [k for k in range(NP) if k not in fixed]

    def f(x):
        th = base.copy()
        th[free] = x
        A = to_A(th)
        val, GA = 0.0, np.zeros((K, K))
        for xb, xa, c2, dt in pairs:
            mu = c2 * expm(dt * A) @ xb
            if not np.all(np.isfinite(mu)) or np.any(mu < -1e-9):
                return 1e300, np.zeros(len(free))
            mu = np.maximum(mu, 1e-300)
            val += float(np.sum(mu - xa * np.log(mu)))
            _, L = expm_frechet(dt * A.T, c2 * (1.0 - xa / mu) @ xb.T)
            GA += dt * L
        return val, np.einsum("pij,ij->p", DA, GA)[free]

    x0 = base[free].copy()
    bnds = [(None, None) if k < K else (0.0, None) for k in free]
    x0 = np.array([max(v, 0.0) if k >= K else v for v, k in zip(x0, free)])
    try:
        r = minimize(f, x0, jac=True, method="L-BFGS-B", bounds=bnds,
                     options={"maxiter": 20000, "ftol": 1e-15, "gtol": 1e-10})
        x = r.x
    except Exception:
        x = x0
    th = base.copy()
    th[free] = x
    return th


def rate_condition(name):
    """Rate model for one condition (runs in its own process).
    Search: closed form, several starts. Values: SIRPY. Every reported
    likelihood is computed from SIRPY solutions and polished with SIRPY."""
    SOLVES["solves"] = SOLVES["verified"] = 0      # count this task only
    data, _ = load_data()
    pairs = data[name]
    fp = Fates([pairs])
    P_hat = fp.split(fp.solve())[0]
    dt = pairs[0][3]
    lik = RateLik(pairs)
    rng = np.random.default_rng(20261011)

    def perturbed(th, fixed=()):
        p = th.copy()
        p[:K] += rng.normal(0, 0.02, K)
        p[K:] = np.maximum(p[K:], 1e-3) * np.exp(rng.normal(0, 0.4, NP - K))
        for k in fixed:
            p[k] = 0.0
        return p

    def cf_value(th):
        val = 0.0
        A = to_A(th)
        for xb, xa, c2, dt_ in pairs:
            mu = np.maximum(c2 * expm(dt_ * A) @ xb, 1e-300)
            val += float(np.sum(mu - xa * np.log(mu)))
        return val

    # ---- full fit: closed-form search from several starts, SIRPY polish
    lin = (P_hat - np.eye(K)) / dt
    s0 = start_from_fates(P_hat, dt)
    s1 = from_A(np.where(np.eye(K, dtype=bool), lin, np.maximum(lin, 1e-3)))
    # forty starts: the likelihood has several local optima (for BMP4 the
    # best one is reached from about one start in five), so the three best
    # candidates are each polished with SIRPY and the best is kept
    starts = [s0, s1] + [perturbed(s0 if i % 2 else s1) for i in range(38)]
    cands = sorted((cf_fit(pairs, s) for s in starts), key=cf_value)
    best = min((lm_fit(lik, c) for c in cands[:3]), key=lambda r: r["nll"])
    th = best["theta"]
    phi = best["pearson"][0] / (best["pearson"][1] - NP)

    # ---- likelihood-ratio test of 'this switch is zero' (one-sided 5%: 2.706)
    # All 30 switches for the untreated network and TMZ; for control and BMP4
    # the switch the paper highlights (state 2 -> state 5).
    scope = LAB[K:] if name in ("untreated", "TMZ") else ["s2->5"]
    nll_k = {}
    for lab in scope:
        k = LAB.index(lab)
        if th[k] <= 1e-10:
            nll_k[lab] = best["nll"]
            continue
        st = th.copy()
        st[k] = 0.0
        tries = [cf_fit(pairs, s, {k: 0.0})
                 for s in [st] + [perturbed(st, (k,)) for _ in range(5)]]
        if 2 * (min(map(cf_value, tries)) - best["nll"]) / phi > 1.5:
            tries += [cf_fit(pairs, perturbed(st, (k,)), {k: 0.0}) for _ in range(6)]
        top = lm_fit(lik, min(tries, key=cf_value), fixed={k: 0.0})   # SIRPY
        nll_k[lab] = top["nll"]
        if top["nll"] < best["nll"] - 1e-6:
            ref = lm_fit(lik, top["theta"])        # better full optimum found
            if ref["nll"] < best["nll"]:
                best, th = ref, ref["theta"]
    lr = {lab: max(2 * (v - best["nll"]) / phi, 0.0) for lab, v in nll_k.items()}
    est = [s for s in scope if lr[s] > 2.706]
    cov = phi * np.linalg.pinv(best["F"], rcond=1e-12)
    Pr = Propagator(to_A(th), dt)                 # closed-form cross-check
    cf = float(np.max(np.abs(Pr.Phi(dt) - expm(dt * to_A(th)))))
    res = {"theta": dict(zip(LAB, th.tolist())), "nll": best["nll"],
           "score": best["score"], "gap": best["gap"], "phi": phi,
           "lr_zero": lr, "tested": list(scope), "established": est,
           "cov": cov.tolist(), "A": to_A(th).tolist(),
           "sirpy_vs_closed_form": cf}
    return name, res, dict(SOLVES)


# ============================================================== main
def main():
    say("SIRPY", sirpy.__version__)
    data, x_day8 = load_data()
    R = {"sirpy_version": sirpy.__version__}

    # ------------------------------------------------ part 1: 8-day fates
    say("Part 1: 8-day fates (no differential equations)")
    joint = Fates([data["control"], data["TMZ"]])
    v_hat = joint.solve()
    phi_j = joint.phi(v_hat)
    Pc, Pt = joint.split(v_hat)

    def e(cond, i, j):
        a = np.zeros(joint.n)
        a[{"control": 0, "TMZ": 1}[cond] * K * K + i * K + j] = 1.0
        return a

    def col(cond, j, rows=range(K)):
        return sum(e(cond, i, j) for i in rows)

    not5 = [0, 1, 2, 3, 5]
    tests = {
        "state2_all_control": col("control", 1),
        "state2_all_TMZ": col("TMZ", 1),
        "state2_all_diff": col("TMZ", 1) - col("control", 1),
        "state2_to5_control": e("control", 4, 1),
        "state2_to5_TMZ": e("TMZ", 4, 1),
        "state2_to5_diff": e("TMZ", 4, 1) - e("control", 4, 1),
        "state2_not5_control": col("control", 1, not5),
        "state2_not5_TMZ": col("TMZ", 1, not5),
        "state2_not5_diff": col("TMZ", 1, not5) - col("control", 1, not5),
        "state1_to5_diff": e("TMZ", 4, 0) - e("control", 4, 0),
        "state4_to5_diff": e("TMZ", 4, 3) - e("control", 4, 3),
        "state5_all_diff": col("TMZ", 4) - col("control", 4),
    }
    fates = {}
    for k, a in tests.items():
        fates[k] = joint.interval(a, v_hat, phi_j)
        say(f"   {k}: {fates[k][0]:+.3f}  95% [{fates[k][1]:.3f}, {fates[k][2]:.3f}]")
    R["fates"] = {"phi": phi_j, "P_control": Pc.tolist(), "P_TMZ": Pt.tolist(),
                  "intervals": fates}

    # ------------------------------------------------ part 2: rates (SIRPY)
    names = ("untreated", "control", "TMZ", "BMP4")
    workers = max(1, min(len(names), os.cpu_count() or 1))
    say(f"Part 2: rates with SIRPY ({workers} conditions at a time; "
        "this is the long part)")
    rates = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name, res, solves in pool.map(rate_condition, names):
            rates[name] = res
            for k in SOLVES:
                SOLVES[k] += solves[k]
            say(f"   {name}: -loglik {res['nll']:.6f}, score {res['score']:.1e}, "
                f"overdispersion {res['phi']:.2f}; switches established "
                f"{len(res['established'])} of {len(res['tested'])} tested: "
                f"{res['established']}")
    R["rates"] = rates

    # ------------------------------------------------ part 3: order (SIRPY)
    say("Part 3: treatment order with SIRPY")
    x0 = x_day8 / x_day8.sum()                 # day-8 population, as fractions
    drugs = {"control": "control", "TMZ": "TMZ", "BMP4": "BMP4"}
    th_d = {d: np.array([rates[d]["theta"][l] for l in LAB]) for d in drugs}
    cov_d = {d: np.array(rates[d]["cov"]) for d in drugs}
    prop = {d: Propagator(to_A(th_d[d]), 16.0) for d in drugs}
    one = np.ones(K)

    def schedule(blocks):
        """blocks: [(drug, days), ...]. Returns fold expansion, state-5 share,
        and the gradient of log(fold) per drug (for the delta method)."""
        ys = [x0]
        for d, t in blocks:
            ys.append(prop[d].Phi(t) @ ys[-1])
        f = float(one @ ys[-1])
        grads = {d: np.zeros(NP) for d in drugs}
        back = one.copy()
        for n in range(len(blocks) - 1, -1, -1):
            d, t = blocks[n]
            if t > 0:
                GA = prop[d].grad_A(np.outer(back, ys[n]), t)
                grads[d] += np.einsum("pij,ij->p", DA, GA) / f
            back = prop[d].Phi(t).T @ back
        return f, ys[-1][4] / f, grads

    def interval(grads, f):
        sd = math.sqrt(sum(g @ cov_d[d] @ g for d, g in grads.items()))
        return [f, f * math.exp(-1.96 * sd), f * math.exp(1.96 * sd)]

    # 16 days built from the 8-day exposures that were actually measured
    plans = {"control, control": [("control", 8), ("control", 8)],
             "TMZ, TMZ": [("TMZ", 8), ("TMZ", 8)],
             "BMP4, BMP4": [("BMP4", 8), ("BMP4", 8)],
             "TMZ, then BMP4": [("TMZ", 8), ("BMP4", 8)],
             "BMP4, then TMZ": [("BMP4", 8), ("TMZ", 8)]}
    order = {"plans": {}}
    for lab, blocks in plans.items():
        f, s5, g = schedule(blocks)
        order["plans"][lab] = {"fold": interval(g, f), "share5": s5}
    fa, _, ga = schedule(plans["TMZ, then BMP4"])
    fb, _, gb = schedule(plans["BMP4, then TMZ"])
    gr = {d: ga[d] - gb[d] for d in drugs}
    order["ratio_8_8"] = interval(gr, fa / fb)
    # other switch days: model extrapolation between the measured 8-day points
    order["switch_day"] = {}
    for tau in (2, 4, 6, 10, 12, 14):
        f, s5, g = schedule([("TMZ", tau), ("BMP4", 16 - tau)])
        order["switch_day"][tau] = interval(g, f)
    order["growth_per_day"] = {n: float(max(np.linalg.eigvals(np.array(rates[n]["A"])).real))
                               for n in ("control", "TMZ", "BMP4", "untreated")}
    R["order"] = order
    r8 = order["ratio_8_8"]
    say(f"   TMZ 8 d then BMP4 8 d, versus the reverse: ratio {r8[0]:.3f} [{r8[1]:.3f}, {r8[2]:.3f}]")

    R["sirpy_solves"] = dict(SOLVES)
    json.dump(R, open(os.path.join(OUT, "results.json"), "w"), indent=1)
    figures(R)
    report(R)


# ============================================================== figures
INK, RED, BLUE, GREY, GREEN, GOLD = "#14304a", "#c0392b", "#2e6da4", "#6b7785", "#2e7d4f", "#e0a800"


def figures(R):
    plt.rcParams.update({"font.size": 12})
    F = R["fates"]["intervals"]
    # ---- figure 1: kill or switch
    fig, ax = plt.subplots(figsize=(7.6, 5.0), dpi=200)
    groups = [("descendants in state 5", "state2_to5"), ("descendants in other states", "state2_not5")]
    for gi, (lab, key) in enumerate(groups):
        for ci, (cond, colr) in enumerate((("control", GREY), ("TMZ", RED))):
            est, lo, hi = F[f"{key}_{cond}"]
            x = gi * 3 + ci
            ax.bar(x, est, color=colr, width=0.8, label=cond if gi == 0 else None)
            ax.errorbar(x, est, yerr=[[est - lo], [hi - est]], color=INK, capsize=6, lw=1.6)
    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel("descendants per state-2 cell, 8 days later")
    ax.set_title("Temozolomide and state-2 cells: kill or switch?\n"
                 "state-5 descendants unchanged, the others lost", fontsize=12.5, color=INK)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.text(0.5, 0.005, "95% intervals. Data: Larsson et al., Mol Syst Biol 2021 (U3065MG, barcoded scRNA-seq).",
             ha="center", fontsize=8.5, color=GREY)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(OUT, "fig1_kill_or_switch.png"))
    plt.close(fig)

    # ---- figure 2: which switches the data establish
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2), dpi=200)
    ang = np.pi / 2 - 2 * np.pi * np.arange(K) / K
    pos = np.c_[np.cos(ang), np.sin(ang)]
    for ax, name, title in ((axes[0], "untreated", "Untreated (the paper's main network)"),
                            (axes[1], "TMZ", "Temozolomide")):
        r = R["rates"][name]
        for i, j in OFF:
            ax.plot(*zip(pos[j], pos[i]), color="#d5dbe1", lw=0.8, zorder=1)
        for s in r["established"]:
            j, i = (int(v) - 1 for v in s[1:].split("->"))
            ax.add_patch(FancyArrowPatch(pos[j] * 0.86, pos[i] * 0.86, arrowstyle="-|>",
                                         mutation_scale=16, lw=2.2, color=GREEN,
                                         connectionstyle="arc3,rad=0.12", zorder=3))
        for k in range(K):
            ax.scatter(*pos[k], s=900, color="white", edgecolor=INK, lw=2, zorder=4)
            ax.text(*pos[k], str(k + 1), ha="center", va="center", fontsize=13,
                    weight="bold", color=INK, zorder=5)
        just = [s for s in r["established"] if r["lr_zero"][s] < 4.0]
        sub = f" ({len(just)} only just)" if just else ""
        ax.set_title(f"{title}\n{len(r['established'])} of 30 switches established{sub}",
                     fontsize=12, color=INK)
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.set_aspect("equal")
        ax.axis("off")
    fig.text(0.5, 0.02, "Green arrow: the data rule out 'this switch is zero' (likelihood-ratio test, 5%). "
             "Grey: possible switches the data do not establish.", ha="center", fontsize=9, color=GREY)
    fig.savefig(os.path.join(OUT, "fig2_switches.png"))
    plt.close(fig)

    # ---- figure 3: treatment order (16 days from the measured 8-day blocks)
    O = R["order"]["plans"]
    order = ["control, control", "BMP4, BMP4", "TMZ, TMZ", "BMP4, then TMZ", "TMZ, then BMP4"]
    labels = {"control, control": "no treatment",
              "BMP4, BMP4": "BMP4 for 16 days",
              "TMZ, TMZ": "TMZ for 16 days",
              "BMP4, then TMZ": "BMP4 8 days, then TMZ 8 days",
              "TMZ, then BMP4": "TMZ 8 days, then BMP4 8 days"}
    cols = {"control, control": GREY, "BMP4, BMP4": BLUE, "TMZ, TMZ": RED,
            "BMP4, then TMZ": "#7d8fb3", "TMZ, then BMP4": GREEN}
    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=200)
    for y, key in enumerate(order):
        f, lo, hi = O[key]["fold"]
        ax.barh(y, f, color=cols[key], height=0.62)
        ax.errorbar(f, y, xerr=[[f - lo], [hi - f]], color=INK, capsize=5, lw=1.5)
        ax.text(hi + 1.0, y, f"{f:.0f}x, state 5: {100 * O[key]['share5']:.0f}%",
                va="center", fontsize=10, color=INK)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([labels[k] for k in order])
    ax.set_xlabel("cells after 16 days, per cell at the start")
    ax.set_xlim(0, max(O[k]["fold"][2] for k in order) * 1.35)
    ax.set_title("Same two drugs, different order", fontsize=12.5, color=INK, loc="left")
    ax.grid(axis="x", alpha=0.25)
    fig.text(0.5, 0.01, "Model predictions with SIRPY from the 8-day experiments; bars: approximate 95% intervals.",
             ha="center", fontsize=8.5, color=GREY)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(OUT, "fig3_order.png"))
    plt.close(fig)


def report(R):
    F = R["fates"]["intervals"]
    say("done; figures and results.json in", OUT)
    p = lambda k: f"{F[k][0]:.2f} (95% {F[k][1]:.2f} to {F[k][2]:.2f})"
    print("\nNUMBERS FOR THE WRITE-UP (from this computer's run)")
    print(f"  SIRPY solves: {R['sirpy_solves']['solves']:,}, verified: {R['sirpy_solves']['verified']:,}")
    print(f"  one state-2 cell, 8 days, all descendants: control {p('state2_all_control')}, TMZ {p('state2_all_TMZ')}")
    print(f"  ... descendants in state 5: control {p('state2_to5_control')}, TMZ {p('state2_to5_TMZ')}, difference {p('state2_to5_diff')}")
    print(f"  ... descendants in other states: difference {p('state2_not5_diff')}")
    print(f"  one state-1 cell, descendants in state 5, TMZ minus control: {p('state1_to5_diff')}")
    print(f"  one state-5 cell, all descendants, TMZ minus control: {p('state5_all_diff')}")
    for n, r in R["rates"].items():
        near = {s: round(v, 2) for s, v in r["lr_zero"].items() if 2.0 < v < 3.5}
        print(f"  {n}: switches established {len(r['established'])} of "
              f"{len(r['tested'])} tested {r['established']}; near the line (LR 2 to 3.5): {near}; "
              f"switch 2->5 LR {r['lr_zero']['s2->5']:.2f} (needs > 2.71); "
              f"score {r['score']:.1e}; SIRPY vs closed form {r['sirpy_vs_closed_form']:.1e}")
    O = R["order"]
    for k, v in O["plans"].items():
        f, lo, hi = v["fold"]
        print(f"  16 days, {k}: {f:.1f} times the cells (95% {lo:.1f} to {hi:.1f}), "
              f"state 5 {100 * v['share5']:.0f}%")
    r = O["ratio_8_8"]
    print(f"  TMZ first versus BMP4 first (8 + 8 days): ratio {r[0]:.3f} (95% {r[1]:.3f} to {r[2]:.3f})")
    sw = {t: f"{v[0]:.1f} ({v[1]:.1f} to {v[2]:.1f})" for t, v in O["switch_day"].items()}
    print(f"  other switch days, TMZ first (model extrapolation): {sw}")
    print("  long-run growth per day: " + ", ".join(f"{k} {v:.3f}" for k, v in O["growth_per_day"].items()))


if __name__ == "__main__":
    main()
