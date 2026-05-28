import pandas as pd


RTS = "/home/ophir/MFLab/rts/output/RTSresult.csv"
MFSIM = "/home/ophir/MFLab/MFSim-cases/dev/validacoes_RTS/3D_Bordbar/output/RTS_results.csv"

rts = pd.read_csv(RTS)#[["i", "j", "k", "x", "y", "z", "T_energy"]]
mfsim = pd.read_csv(MFSIM)#[["i", "j", "k", "x", "y", "z", "T_energy"]]
comp = pd.merge(rts, mfsim, on = ["i", "j", "k", "x", "y", "z"], how = "outer")

def filter_error(df, field, error = 1e-9):
    df[f"error_{field}"] = (df[f"{field}_x"] - df[f"{field}_y"]).abs() / df[f"{field}_x"]
    return df[((df[f"{field}_x"] - df[f"{field}_y"]).abs()) / df[f"{field}_x"] > error]

diff_T = filter_error(comp, "T_energy", 1e-9)
diff_Srad = filter_error(comp, "S_rad", 1e-9)
diff_Grad = filter_error(comp, "G", 1e-9)
diff_Qradw = filter_error(comp, "q_radw", 1e-9)