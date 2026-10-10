# Theophylline after one oral dose — two answers, one curve

A worked example of SIRPY fitting a classic pharmacokinetic model and reporting what the data can and cannot decide.

## Source

Data from the R dataset `Theoph`. Boeckmann, Sheiner and Beal (1994) report a study by Dr. Robert Upton of the anti-asthma drug theophylline. We use subject 1: one oral dose of 4.02 mg/kg, then 11 serum concentrations over 24.4 hours. All 11 samples are used as recorded, including the pre-dose sample (0.74 mg/L at t = 0).

## Data

``` python
import numpy as np
t = np.array([0, 0.25, 0.57, 1.12, 2.02, 3.82, 5.1, 7.03, 9.05, 12.12, 24.37])    # hours
c = np.array([0.74, 2.84, 6.57, 10.5, 9.66, 8.58, 8.36, 7.47, 6.89, 5.94, 3.28])  # mg/L
dose = 4.02                                                                        # mg/kg
```

In R: `Theoph[Theoph$Subject == 1, ]`

## Model

The standard one-compartment model with first-order absorption. $G(t)$ is the drug not yet absorbed (mg/kg) and $C(t)$ is the concentration in serum (mg/L):

$$\frac{dG}{dt} = -k_a G, \qquad G(0) = 4.02$$

$$\frac{dC}{dt} = \frac{k_a}{V}G - k_e C, \qquad C(0) = 0, \qquad V = \frac{CL}{k_e}$$

The unknowns are the absorption rate $k_a$, the elimination rate $k_e$ and the clearance $CL$, estimated on the log scale by maximum likelihood with normal errors of constant variance (the same optimum as least squares). As usual for oral data, $CL$ and $V$ are apparent values (divided by the unknown bioavailability).

## Result

**Two answers fit equally well.** Both have the same sum of squared errors, 4.28600902 (residual standard error 0.732 mg/L), and their curves agree to 10 decimal places.

| | Answer A: absorbed fast, cleared slowly | Answer B: absorbed slowly, cleared fast |
| :--- | :--- | :--- |
| absorption rate ka (1/h) | 1.777 [1.351, 2.362] | 0.0540 [0.0484, 0.0609] |
| elimination rate ke (1/h) | 0.0540 [0.0484, 0.0609] | 1.777 [1.351, 2.362] |
| clearance CL (L/h per kg) | 0.0199 [0.0181, 0.0219] | 0.0199 [0.0181, 0.0219] |
| volume V (L/kg) | 0.369 | 0.0112 |
| absorption half-life | 23 min | 12.8 h |
| elimination half-life | 12.8 h | 23 min |

Brackets: 95% profile-likelihood intervals from SIRPY's estimator.

**What the data decide.** Clearance: 0.0199 L/h per kg, 95% interval [0.0181, 0.0219], the same in both answers.

**What the data cannot decide.** Which rate is absorption and which is elimination. The textbook solution, with dose $D$ = 4.02 mg/kg, shows why:

$$C(t) = \frac{D k_a k_e}{CL (k_a - k_e)}\left(e^{-k_e t} - e^{-k_a t}\right)$$

It does not change when $k_a$ and $k_e$ are swapped with $CL$ held fixed, so no fit to these samples can tell the two answers apart. Pharmacologists call this flip-flop.

**What settles it.** Answer B needs $V$ = 0.0112 L/kg, far below the volume of the blood plasma alone (about 0.04 L/kg). The drug is measured in serum, so its volume of distribution cannot be smaller than the plasma volume, and a bioavailability below 100% would only make B's true volume smaller. Answer B is impossible; answer A is the real one.

**Why it matters.** For this formulation both answers predict the same concentrations. For an IV dose or a new formulation they predict very different ones.

## How the analysis was done

1. **Map.** A 29 × 29 grid of $(k_a, k_e)$, evenly spaced on a log scale from 0.01 to 10 per hour, one SIRPY solve per pair with the best clearance for that pair (841 solves). The misfit has two equally deep valleys, mirror images across the line $k_a = k_e$ (Figure 2).
2. **Fit.** SIRPY's estimator, started from the best grid point in each valley, with 95% profile-likelihood intervals (18,438 solves for A and 18,305 for B).

## How we know the numbers are right

- 37,584 ODE solves, every one returned `solved+verified` (a failed solve is discarded, never used).
- The exact textbook formula above agrees with SIRPY's curves to 11 decimal places at the sample times.
- R's standard fit of the same data matches every printed digit (next section).

## Check it yourself in R

Install R (free, from r-project.org) and run

``` r
example(SSfol)
```

R fits the same model to the same 11 samples and prints lKe = −2.9196, lKa = 0.5752 and lCl = −3.9159 (natural logs of $k_e$, $k_a$ and $CL$), with residual standard error 0.732 on 8 degrees of freedom. These are answer A: $e^{-2.9196} = 0.0540$, $e^{0.5752} = 1.777$, $e^{-3.9159} = 0.0199$.

## Exact call

``` python
import sirpy

r = sirpy.solve_ivp(
    system={"Dy1": "-ka*y1",             # G: drug not yet absorbed
            "Dy2": "ka*y1/V - ke*y2"},   # C: concentration in serum
    x0=0, x1=24.37, ic={"y1": 4.02, "y2": 0.0},
    params={"ka": 1.77741, "ke": 0.053955, "V": 0.36926},   # answer A
    verbose=False)

r["status"]                  # 'solved+verified'
r["evaluate"]("y2", 2.02)    # concentration at 2.02 hours
```

For answer B, use `params={"ka": 0.053955, "ke": 1.77741, "V": 0.011209}`.



## Reproducibility

Every SIRPY number and figure above comes from runs on our own computers, and the R check can be run by anyone. The solver engine is proprietary and not distributed here.

## Credit

`Theoph.csv` is the R dataset Theoph (Boeckmann, Sheiner and Beal, 1994), included for convenience. It is third-party data, not covered by LQI’s copyright.

Data: Boeckmann, A. J., Sheiner, L. B. and Beal, S. L. (1994). *NONMEM Users Guide: Part V.* NONMEM Project Group, University of California, San Francisco. Distributed with R as the dataset `Theoph`.