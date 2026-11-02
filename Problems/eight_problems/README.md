# Eight problems, pure SIRPY

Each problem below is stated in full, followed by the one call a user types and what SIRPY printed back. The quoted output is SIRPY's own (version 0.67.0), unedited except where a cut is marked `[…]`; the complete console output is in [`sirpy_pure_output.txt`](sirpy_pure_output.txt). The script that produces it is [`sirpy_pure.py`](sirpy_pure.py); it uses only the public API and writes nothing to disk. A second script, [`sirpy_showcase.py`](sirpy_showcase.py), adds what SIRPY does not know: the exact solutions for comparison and a scipy run next to problem 2.

Two of SIRPY's lines appear on most problems:

- **VERIFIED** is the residual check: the delivered function is substituted back into the equation at points away from the piece boundaries, and the leftover is measured.
- **ERROR BOUND** is a proven statement about the distance to the exact solution: existence, uniqueness, and the bound. When the proof does not close, SIRPY says so and says where.

---

## 1. y′ = sin(y), y(0) = 1, on [0, 3]

```python
r = sirpy.solve_ivp(order=1, f_str="sin(y)", x0=0, x1=3, ic={0: 1.0})
print(sirpy.error_bound_text(r.error_bound))
```

```
INFO: run summary — series degree 20, iterations=20 (sirpy default), 6 series segment(s).
INFO: verification by substitution — VERIFIED: satisfies the equation to 1.27e-11 on [0, 3]. Checked by exact coefficient differentiation at interior points, not at the interface nodes. Interface continuity across 5 node(s): max one-sided jump in (y..y0) = 0.0e+00 — each side's own coefficients evaluated exactly at the node. Also checked: initial data to 0.00e+00.
ERROR BOUND (proven): the problem has exactly one solution on [0, 3] and the delivered function differs from it by at most 8.761e-11 in every component (Euclidean norm); relative to the solution scale 2.96 that is 2.96e-11; largest proven defect on a piece 3.74e-10; 6 piece(s).
```

The exact solution is 2 arctan(tan(1/2) eˣ); the delivered function is within 8.9e-12 of it (`sirpy_showcase.py`), inside the proven 8.761e-11.

## 2. y′ = y², y(0) = 1, asked to x = 2

The solution 1/(1 − x) blows up at x = 1, inside the requested interval.

```python
r = sirpy.solve_ivp(order=1, f_str="y**2", x0=0, x1=2, ic={0: 1})
print(sirpy.error_bound_text(r.error_bound))
```

```
INFO: verification by substitution — VERIFIED: satisfies the equation to 3.62e-02 on [0, 0.999944] (the equation's scale reaches 3.14e+08; max RELATIVE residual, taken pointwise, is 1.15e-10). […]
WARNING: the solution appears to blow up near x ~ 0.999944 (the continuation segments shrank to a point — the signature of a pole/finite-time blow-up). The solution does not exist past that point; it is returned on [0, 0.999944] over 30 segments, and values requested beyond are not meaningful.
INFO: a singularity of the solution is estimated near x = 1 (nearest pole; two independent estimators agree). The series is trustworthy only for |x - 0| < 1.
WARNING: your interval reaches x = 2, PAST the estimated singularity at 1. Values beyond it are not the same solution branch; the result is meaningful only up to the singularity.
ERROR BOUND (proven): the problem has exactly one solution on [0, 0.999944] and the delivered function differs from it by at most 2.504e-03 in every component (max norm); relative to the solution scale 1.77e+04 that is 1.41e-07; largest proven defect on a piece 3.62e-02; 30 piece(s).
```

The absolute bound is large because the solution is large there (1.77e+04); relative to the solution scale it is 1.4e-07. For comparison, scipy's RK45 on the same problem stops at x = 0.999929 with "Required step size is less than spacing between numbers" (`sirpy_showcase.py`).

## 3. Torricelli's draining tank, h′ = −k√h, h(0) = 11.1 cm, k = 0.0025812·√981, to t = 100 s

The tank empties at t = 82.42 s. From then on the equation is not Lipschitz at h = 0 and has infinitely many solutions; the physical one is h ≡ 0.

```python
r = sirpy.solve_ivp(order=1, f_str="-0.0025812*sqrt(981)*sqrt(y)", x0=0, x1=100, ic={0: 11.1})
print(sirpy.error_bound_text(r.error_bound))
```

```
INFO: verification by substitution — VERIFIED: satisfies the equation to 8.18e-12 on [0, 100]. […]
INFO: BRANCH POINT of the equation at x = 82.42050169: the square root's radicand vanishes at y = 0, where the equation is not Lipschitz. The classical maximal solution continues as the CONSTANT y = 0 up to x = 100, and that glued solution is what is returned and verified.
ERROR BOUND: not closed -- the tube argument did not close on piece 28 starting at x = 82.413: the logarithmic norm of the Jacobian is too large there for a Lipschitz-type bound (strong stiffness, or no Lipschitz constant at all -- e.g. 1/sqrt(y) at a square-root branch point). On [0, 82.413) the solution is unique and within 6.134e-03 of the delivered function (28 of 30 pieces); beyond that VERIFIED means the residual check. The proven bound is below 1e-6 of the solution scale on [0, 82.2253]. Assumptions: a constant that is not a rational or a binary64 number was evaluated by SymPy to 50 significant digits and widened by 1e-45 relative (not an interval-arithmetic proof).
```

What "verified" means here, part by part: on [0, 82.413) a residual check and a proven bound; on [82.413, 100] the delivered solution is the constant branch h = 0, verified by residual check only (h′ = 0 = −k√0). No Lipschitz bound exists across the branch point, so none is claimed. The closed form gives the emptying time 2√11.1/k = 82.420515 s (printed by `sirpy_showcase.py`); SIRPY located it from the equation alone as 82.42050169.

## 4. Van der Pol, y″ = 1000 (1 − y²) y′ − y, y(0) = 2, y′(0) = 0, on [0, 10]

A stiff problem. No Jacobian is supplied and no tolerance is tuned.

```python
r = sirpy.solve_ivp(order=2, f_str="1000*(1-y**2)*y1 - y", x0=0, x1=10, ic={0: 2, 1: 0})
print(sirpy.error_bound_text(sirpy.error_bound(r)))   # explicit call: above 2000 pieces the bound is on request
```

```
INFO: run summary — series degree 20, iterations=20 (sirpy default), 3388 series segment(s).
INFO: verification by substitution — VERIFIED: satisfies the equation to 2.09e-07 on [0, 10]. Checked by exact coefficient differentiation at interior points, not at the interface nodes. Interface continuity across 3387 node(s): max one-sided jump in (y..y1) = 4.1e-11 — each side's own coefficients evaluated exactly at the node. Also checked: initial data to 0.00e+00.
[three run-management notices omitted; see sirpy_pure_output.txt]
ERROR BOUND (proven): the problem has exactly one solution on [0, 10] and the delivered function differs from it by at most 3.219e-05 in every component (max norm, weights varying along the solution, floor 0.01); relative to the solution scale 2 that is 1.61e-05; largest proven defect on a piece 1.24e-07; 3388 piece(s).
```

## 5. Bratu, y″ = −eʸ, y(0) = 0, y(1) = 0 — two solutions

This boundary value problem has two solutions. Each starting slope reaches one of them; both are verified.

```python
r = sirpy.solve_bvp(order=2, f_str="-exp(y)", x0=0, x1=1, left_bc={0: 0, 1: "gamma"}, right_bc={0: 0}, gamma0=0.5)
r = sirpy.solve_bvp(order=2, f_str="-exp(y)", x0=0, x1=1, left_bc={0: 0, 1: "gamma"}, right_bc={0: 0}, gamma0=11.0)
```

```
INFO: verification by substitution — VERIFIED: satisfies the equation to 7.93e-11 on [0, 1], and the boundary data to 1.37e-13. […]
INFO: exact first-integral test — the motion is bounded (a turning point is reached); this solution has no blow-up.
INFO: verification by substitution — VERIFIED: satisfies the equation to 2.31e-09 on [0, 1] (the equation's scale reaches 59.8; max RELATIVE residual, taken pointwise, is 6.12e-11), and the boundary data to 6.12e-15. […]
INFO: exact first-integral test — the motion is bounded (a turning point is reached); this solution has no blow-up.
```

The two slopes found are y′(0) = 0.5493527288 and y′(0) = 10.84689902, with y(1/2) = 0.14053921 and 4.0914672 (`sirpy_showcase.py`).

## 5b. y″ = −y, y(0) = 1, y(π) = 1 — no solution exists

```python
r = sirpy.solve_bvp(order=2, f_str="-y", x0=0, x1=math.pi, left_bc={0: 1, 1: "gamma"}, right_bc={0: 1})
```

```
CERTIFICATE OF NONEXISTENCE -- no slope satisfies these boundary conditions.
INFO: A scan of 35 trial slopes gamma in [-32, 32] gives a terminal mismatch y(3.141592653589793) - (1) that is NEGATIVE for every trial (trajectories that blow up depart on the same side): the mismatch never changes sign, so no root exists in the scanned range.
INFO: Smallest achievable |mismatch| = 2.000e+00 at gamma = -9.929163e+09; that floor trajectory is itself a verified solution of the ODE (max relative defect 6.3e-07) -- it simply cannot meet y(3.141592653589793) = 1.
INFO: Solutions of this two-point problem exist only up to x1 = b* ~ 3.079113: a saddle-node fold where two solution branches merge (slope gamma* ~ 32.00000). Your x1 = 3.141592653589793 lies 0.0625 beyond the fold. Any x1 <= b* is solvable, with TWO branches; seed gamma0= to select one.
[…]
INFO: Certificate scope: single shooting unknown, scanned range as above (about 4.5x the natural slope scale); a root outside that range would evade the scan. This is a measured finding, not a raised error: inspect r['certificate'] for the numbers.
```

Every solution of y″ = −y is y = cos x + c sin x, so y(π) = −1 for every c: no solution exists, as certified.

## 6. Volterra integral equation, y(t) = sin(2t)/2 + (3/2) ∫₀ᵗ y(s)² cos(t − s) ds, on [0, 3]

```python
r = sirpy.solve_volterra(phi="sin(2*t)/2", kernel="cos(t-s)", f="y**2", lam=1.5, a=0, t1=3)
print(sirpy.error_bound_text(r.error_bound))
```

```
INFO: the equation supplies its own initial data: y(0) = phi(0) = 0,  y'(0) = phi'(0) + lam*K(0,0)*f = 1
INFO: verification by substitution — VERIFIED: the returned series satisfies the Volterra equation to 1.19e-13 (relative 1.03e-13). The right-hand side was recomputed by quadrature from the series itself. […]
ERROR BOUND (proven): the problem has exactly one solution on [0, 3] and the delivered function differs from it by at most 5.241e-10 in every component (Euclidean norm, weights varying along the solution, floor 0.01); relative to the solution scale 1.33 that is 3.93e-10; largest proven defect on a piece 5.70e-12; 5 piece(s).
```

The exact solution is sin t; the delivered function is within 4.7e-13 of it, inside the proven 5.241e-10 (`sirpy_showcase.py`).

## 7. Parameter estimation, Ornstein–Uhlenbeck dX = κ(m − X) dt + σ dW

1,501 weekly observations simulated with κ = 2, m = 0.05, σ = 0.02 (the generating code is in `sirpy_pure.py`).

```python
r = sirpy.fit_mle("ou", np.array(x), dt)
```

```
INFO: MLE by certified score root-finding — 9 start(s) across the box converged to 1 mode(s); agreement across starts is UNIQUENESS EVIDENCE within the box, not a proof.
INFO: verification — VERIFIED: max |score| = 0.00e+00 at theta-hat; observed information is negative-definite Hessian (all curvature directions checked); all parameters identified.
INFO: theta-hat with SE and 95% profile-likelihood CI — kappa = 2.08692 (SE 0.38, CI [1.3394, 2.84535]); m = 0.0485885 (SE 0.0018, CI [0.0448286, 0.0523921]); sigma = 0.0202129 (SE 0.00038, CI [0.0194975, 0.0209739]). CIs by likelihood-ratio bisection, not Wald-only.
INFO: coverage — CHECKED: the score at theta-hat, Hessian curvature in every direction, identifiability (flat combinations named), multi-start mode count, profile-likelihood CIs. NOT checked: model adequacy (goodness-of-fit, residual diagnostics, information criteria, out-of-sample) — a fitted model can be verified as a MAXIMUM and still be the wrong model for the data; the workflow tier for that is on the roadmap.
```

## 8. The client document

`sirpy.certificate_text(r)` renders the verification report for any result; `sirpy_pure.py` prints it for problem 1. Section 2 states the finding and, when a proof exists, the error-bound theorem in full.

---

## Running it

From the folder that holds the `sirpy` package:

```
python sirpy_pure.py        # the eight problems, SIRPY's own output only (about a minute)
python sirpy_showcase.py    # the same problems with exact solutions and comparisons added
```

SIRPY is offered as a service by Luminary Quantum Institute LLC. The scripts above use only the public API; the engine is not distributed.
