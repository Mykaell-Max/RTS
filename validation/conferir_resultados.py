import pandas as pd
from pathlib import Path
import numpy as np

ROOT = Path("/home/ophir/MFLab/validacoes_RTS")

def is_edge_3D(i,j,k, nxi, nyi, nzi):
    i = i == 1 or i == nxi
    j = j == 1 or j == nyi
    k = k == 1 or k == nzi

    return (i and j) or (i and k) or (j and k)

def is_edge_2D(i,j, nxi, nyi):
    i = i == 1 or i == nxi
    j = j == 1 or j == nyi

    return (i and j)

def is_edge_1D(i, nxi):
    i = i == 1 or i == nxi

    return i

def compare_csv_3D(case_name: str, nx, ny, nz):
    case_root = Path(ROOT, case_name)
    rts_path = Path(case_root, "RTS", "output", "RTSresult.csv")
    mfsim_path = Path(case_root, "MFSim", "output", "RTS_results.csv")

    rts = pd.read_csv(rts_path)
    mfsim = pd.read_csv(mfsim_path)

    rts = rts[rts.apply(lambda r: not is_edge_3D(r.i, r.j, r.k, nx, ny, nz), axis = 1)]
    mfsim = mfsim[mfsim.apply(lambda r: not is_edge_3D(r.i, r.j, r.k, nx, ny, nz), axis = 1)]

    pd.testing.assert_frame_equal(
        rts, mfsim, rtol = 1e-5, atol = 1e-8, check_dtype = False
    )

    print("APROVADO")

def compare_csv_2D(case_name: str, nx, ny):
    case_root = Path(ROOT, case_name)
    rts_path = Path(case_root, "RTS", "output", "RTSresult.csv")
    mfsim_path = Path(case_root, "MFSim", "output", "RTS_results.csv")

    rts = pd.read_csv(rts_path)
    mfsim = pd.read_csv(mfsim_path)

    rts = rts[rts.k == 2]
    mfsim = mfsim[mfsim.k == 2]

    rts = rts[rts.apply(lambda r: not is_edge_2D(r.i, r.j, nx, ny), axis = 1)]
    mfsim = mfsim[mfsim.apply(lambda r: not is_edge_2D(r.i, r.j, nx, ny), axis = 1)]

    pd.testing.assert_frame_equal(
        rts, mfsim, rtol = 1e-5, atol = 1e-8, check_dtype = False
    )

    print("APROVADO")

def compare_csv_1D(case_name: str, nx):
    case_root = Path(ROOT, case_name)
    rts_path = Path(case_root, "RTS", "output", "RTSresult.csv")
    mfsim_path = Path(case_root, "MFSim", "output", "RTS_results.csv")

    rts = pd.read_csv(rts_path)
    mfsim = pd.read_csv(mfsim_path)

    rts = rts[rts.k == 2]
    mfsim = mfsim[mfsim.k == 2]

    rts = rts[rts.j == 2]
    mfsim = mfsim[mfsim.j == 2]

    # rts = rts[rts.apply(lambda r: not is_edge_1D(r.i, nx), axis = 1)]
    # mfsim = mfsim[mfsim.apply(lambda r: not is_edge_1D(r.i, nx), axis = 1)]

    pd.testing.assert_frame_equal(
        rts, mfsim, rtol = 1e-12, atol = 1e-8, check_dtype = False
    )

    print("APROVADO")

#%%
# CASE = ("3D_Bordbar", 19, 36, 19)
# CASE = ("3D_Hsu", 23, 23, 23)
# CASE = ("3D_Soucasse", 44, 44, 44)
# CASE = ("2D_Goutiere", 23, 23)
# CASE = ("2D_Shah", 23, 23, 3)
# CASE = ("2D_Kim", 27, 27)
CASE = ("1D_Bordbar", 30)
# compare_csv_3D(CASE[0], CASE[1], CASE[2], CASE[3])
# compare_csv_2D(CASE[0], CASE[1], CASE[2])
compare_csv_1D(CASE[0], CASE[1])

#%%
case_root = Path(ROOT, CASE[0])
rts_path = Path(case_root, "RTS", "output", "RTSresult.csv")
mfsim_path = Path(case_root, "MFSim", "output", "RTS_results.csv")

rts = pd.read_csv(rts_path)
mfsim = pd.read_csv(mfsim_path)

# VAR = "T_energy"
# VAR = "cappa"
# VAR = "sigma_rad"
# VAR = "beta"
# VAR = "S_rad"
VAR = "G"
# VAR = "XH2O"
# var = "XCO2"

# comp = pd.merge(rts, mfsim, how = "outer", on = ["i", "j", "k"])
# comp = comp[comp.apply(lambda r: not is_edge_2D(r.i, r.j, CASE[1], CASE[2]), axis = 1)] 
comp = comp[comp.apply(lambda r: not is_edge_1D(r.i, CASE[1]), axis = 1)] 
comp = comp[comp.k == 2]
comp = comp[comp.j == 2]
comp["diff%"] = np.abs((comp[f"{VAR}_x"] - comp[f"{VAR}_y"]) / comp[f"{VAR}_x"]) * 100
compf = comp[["i", "j", "k", f"{VAR}_x", f"{VAR}_y"]]
diffs = comp[comp[f"{VAR}_x"] != comp[f"{VAR}_y"]]
diffs
# %%
