# Troesch's Problem — Sample Deliverables

Complete, unedited deliverables from two SIRPY engagements, published so they can be checked independently.

**Equation**

```
y'' = mu*sinh(mu*y),    y(0) = 0,   y(1) = 1
```

A plasma-confinement model, and a standard benchmark: as `mu` grows the solution stays near zero across almost the whole interval, then climbs to 1 inside a boundary layer that thins exponentially.

| Engagement | Finding | Headline |
|---|---|---|
| **mu = 15** | solved + verified | max relative residual **2.98e-10** (criterion 1e-6); boundary data 5.85e-10 |
| **mu = 20** | solved, **NOT verified** | defect confined to a boundary layer ~1.6e-14 wide at x = 1; interior residual 3.6e-15 |

The `mu = 20` engagement is published deliberately. A NOT VERIFIED finding is a measurement, not a fulfillment state — it names where the problem breaks and which lever to pull, and it is delivered in full.

---

## What's here

Each engagement appears at all three tiers (`T1`, `T2`, `T3`), as a separate illustrative order.

| File | What it is |
|---|---|
| `*_artifact.json` | The solution: piecewise power-series coefficients, nodes, verification block, diagnostics, run transcript |
| `*_artifact_json.sha256` | SHA-256 of the artifact |
| `*_ORDER_AND_DELIVERY.txt` | Order manifest: tier, deliverables checklist, artifact binding, reproduction scope, finding |
| `LQI2026Troesch*T1.pdf` | Tier 1 — Verified: report and figures |
| `LQI2026Troesch*T2.pdf` | Tier 2 — Certified Standard: formal certificate, four measurable claims |
| `LQI2026Troesch*T3.pdf` | Tier 3 — Certified Premium: full certificate with diagnostics, localization, claims ledger |
| `LQI2026Troesch*T3_addendum.pdf` | Premium analytics addendum |

Client names and order identifiers are illustrative.

---

## Check the hash

```bash
sha256sum LQI-2026-Troesch-15-T1_artifact.json
```

Expected for Tier 1, mu = 15:

```
60a7f2816691f97aeade4f2b82aa76d33fec944a330531c23c02bbe57e026df3
```

All published digests:

| Artifact | SHA-256 |
|---|---|
| `Troesch-15-T1` | `60a7f2816691f97aeade4f2b82aa76d33fec944a330531c23c02bbe57e026df3` |
| `Troesch-15-T2` | `c751490cd3108fa301eaaffbbbaff6c4bc5b78494a0d801d90f0a1853cc35704` |
| `Troesch-15-T3` | `ce472044612b41fc56e71cb744ee6f1e92d99cd9303a4b79c1f3dfbb1ac3853c` |
| `Troesch-20-T1` | `3e50acb428d92f151d01d3d7f5a1c63887d32551364f55b3bd4c4d2aa772e726` |
| `Troesch-20-T2` | `7bfb7cf8fffc1d8e97d825d1aa30b2a2d30a3b2b57f5531db1cee4516acebc25` |
| `Troesch-20-T3` | `733d06d8c3e373187673aeeb6e7acbaf0f65e98fb79565cfdbc382c4f3b8bd35` |

Each order is a separate engagement, so each artifact carries its own order identifier and therefore its own hash.

A `.gitattributes` guard keeps these files byte-stable across Windows, macOS and Linux checkouts. **If a digest you compute does not match, please tell us** — some download paths alter text files in transit, and we would rather send you the bytes directly.

---

## Rebuild the solution yourself

The artifact carries its own reconstruction rule. On each segment `i`:

```
y(x) = sum_k c_k * (x - segment_centers[i])^k     for x in [nodes[i], nodes[i+1]]
```

with `c_k` from `series_coefficients`. No SIRPY installation is needed — the coefficients, nodes, and the stated equation are sufficient to recompute the reported residuals.

```python
import json
d = json.load(open("LQI-2026-Troesch-15-T1_artifact.json"))
print(d["statement"])
print(d["verification"]["verdict"])
print(d["reconstruction"])
```

---

## How to read the residuals

For `mu = 15` the absolute residual is `1.59e-03` — which looks poor until you see that `sinh(15y)` drives the equation's own terms to `2.45e+07`. Measured against the quantities actually being balanced, the pointwise relative residual is `2.98e-10`.

Every report prints all three: the defect, the scale it must be judged against, and the resulting relative figure. The acceptance criterion is stated before the verdict, not chosen after it.

---

## Scope

We certify that the delivered solution satisfies the **stated** equation, data and constraints. We do not certify that the equation models any particular physical system, and we do not grant regulatory approval. Tier 2 is evidence for a compliance file, not compliance itself.

Not covered at any tier: whether the equation models your system; agreement with another solver (a second opinion, not a reproduction); a solver's tolerance presented as an accuracy claim.

---

## Contact

Pricing depends on the problem. Have a differential equation you think is hard? We would like to see it — we will run it and report honestly what happened.

**info@mitrainstituteofeducation.org** · **https://mitrainstituteofeducation.org/sirpy/**

Please only send problems you are free to share publicly.

---

*SIRPY 0.59.15 · © 2026 Luminary Quantum Institute, LLC (dba Mitra Institute of Education). All rights reserved. The SIRPY solver engine is proprietary and is not distributed here.*
