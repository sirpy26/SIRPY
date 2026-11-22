#!/usr/bin/env python3
"""Optional: rebuild the two CSV files from the authors' original data file.

    pip install mat-io
    python3 extract_data.py expressionMatrices_20210622.mat

The .mat file is in the authors' repository github.com/idalarsson/STAG
(Larsson et al., Mol Syst Biol 2021). It holds MATLAB tables with one row per
cell: cell id, 18-nt lineage barcode, day, treatment, and the authors' state
(1-6). This script writes them unchanged as CSV.
"""
import sys

import numpy as np
import pandas as pd
import matio


def clean(df, statecol):
    return pd.DataFrame({
        "cell": [str(np.ravel(v)[0]) for v in df["Var1"]],
        "barcode": [str(np.ravel(v)[0]) for v in df["X18nt_barcode_sequence"]],
        "day": df["days"].astype(float).astype(int),
        "treatment": df["treatment"].astype(float).astype(int),
        "state": df[statecol].astype(float).astype(int)})


d = matio.load_from_mat(sys.argv[1])
clean(d["expMatOriginal"], "states").to_csv("U3065MG_timecourse_cells.csv", index=False)
clean(d["expMatTreatment"], "StateUpdated").to_csv("U3065MG_treatment_cells.csv", index=False)
print("wrote U3065MG_timecourse_cells.csv and U3065MG_treatment_cells.csv")
