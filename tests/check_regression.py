#!/usr/bin/env python3
"""
check_regression.py  --  compara o output de uma rodada contra o baseline salvo.

Uso:
    python3 tests/check_regression.py <CASO> <novo_csv>

Exemplos:
    python3 tests/check_regression.py 1D_Bordbar validation/1D_Bordbar/RTS/output/RTSresult.csv
    python3 tests/check_regression.py 3D_Hsu     validation/3D_Hsu/RTS/output/RTSresult.csv

Retorna:
    exit 0  se APROVADO (diferença dentro da tolerância)
    exit 1  se REPROVADO (diferença acima da tolerância ou erro de leitura)
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Dimensões de cada caso (usadas para filtrar células de borda)
# ---------------------------------------------------------------------------
CASOS = {
    "1D_Bordbar":  {"dim": 1, "nx": 30},
    "2D_Goutiere": {"dim": 2, "nx": 23, "ny": 23},
    "2D_Kim":      {"dim": 2, "nx": 27, "ny": 27},
    "2D_Shah":     {"dim": 2, "nx": 23, "ny": 23},
    "3D_Bordbar":  {"dim": 3, "nx": 19, "ny": 36, "nz": 19},
    "3D_Hsu":      {"dim": 3, "nx": 23, "ny": 23, "nz": 23},
    "3D_Soucasse": {"dim": 3, "nx": 44, "ny": 44, "nz": 44},
}

# ---------------------------------------------------------------------------
# Funções de filtragem de borda (mesma lógica do conferir_resultados.py)
# ---------------------------------------------------------------------------

def is_edge_3D(i, j, k, nx, ny, nz):
    nxi, nyi, nzi = nx + 2, ny + 2, nz + 2
    ei = (i == 1 or i == nxi)
    ej = (j == 1 or j == nyi)
    ek = (k == 1 or k == nzi)
    return (ei and ej) or (ei and ek) or (ej and ek)

def is_edge_2D(i, j, nx, ny):
    nxi, nyi = nx + 2, ny + 2
    return (i == 1 or i == nxi) and (j == 1 or j == nyi)

# ---------------------------------------------------------------------------
# Comparação principal
# ---------------------------------------------------------------------------

def compare(case_name: str, new_csv: Path, baseline_csv: Path) -> bool:
    if case_name not in CASOS:
        print(f"[ERRO] Caso desconhecido: '{case_name}'")
        print(f"       Casos válidos: {list(CASOS.keys())}")
        return False

    info = CASOS[case_name]
    dim  = info["dim"]

    try:
        novo = pd.read_csv(new_csv)
        ref  = pd.read_csv(baseline_csv)
    except Exception as e:
        print(f"[ERRO] Falha ao ler CSV: {e}")
        return False

    # --- filtrar células internas conforme dimensão ---
    if dim == 3:
        nx, ny, nz = info["nx"], info["ny"], info["nz"]
        novo = novo[~novo.apply(lambda r: is_edge_3D(r.i, r.j, r.k, nx, ny, nz), axis=1)]
        ref  = ref [~ref .apply(lambda r: is_edge_3D(r.i, r.j, r.k, nx, ny, nz), axis=1)]

    elif dim == 2:
        nx, ny = info["nx"], info["ny"]
        novo = novo[novo.k == 2]
        ref  = ref [ref .k == 2]
        novo = novo[~novo.apply(lambda r: is_edge_2D(r.i, r.j, nx, ny), axis=1)]
        ref  = ref [~ref .apply(lambda r: is_edge_2D(r.i, r.j, nx, ny), axis=1)]

    else:  # 1D
        novo = novo[(novo.k == 2) & (novo.j == 2)]
        ref  = ref [(ref .k == 2) & (ref .j == 2)]

    # --- resetar índices para comparação alinhada ---
    novo = novo.reset_index(drop=True)
    ref  = ref .reset_index(drop=True)

    # --- checar mesmo número de linhas ---
    if len(novo) != len(ref):
        print(f"[ERRO] Número de linhas difere: novo={len(novo)} baseline={len(ref)}")
        return False

    if list(novo.columns) != list(ref.columns):
        print(f"[ERRO] Colunas diferem:\n  novo={list(novo.columns)}\n  ref ={list(ref.columns)}")
        return False

    # --- comparação numérica robusta, coluna a coluna ---
    #
    # Observação: o csv_save do RTS usa formato Fortran G30.20. Para números
    # subnormais com expoente de 3 dígitos (ex.: 1e-305) o Fortran imprime
    # "0.123...-305" SEM o "E", o que o pandas lê como string (dtype object).
    # Esses valores sao poeira de ponto flutuante (~0) em células fantasma.
    # Por isso convertemos cada coluna para numérico com errors='coerce'
    # (subnormais malformados viram NaN) e tratamos NaN como 0 antes de
    # comparar com tolerância. Isso é mais robusto que assert_frame_equal,
    # que compara colunas object como texto exato.
    rtol, atol = 1e-5, 1e-8
    falhas = []
    for col in novo.columns:
        a = pd.to_numeric(novo[col], errors="coerce").to_numpy(dtype=float)
        b = pd.to_numeric(ref [col], errors="coerce").to_numpy(dtype=float)
        a = np.nan_to_num(a, nan=0.0)
        b = np.nan_to_num(b, nan=0.0)
        diff = np.abs(a - b)
        tol  = atol + rtol * np.abs(b)
        ruins = diff > tol
        if ruins.any():
            idx = int(np.argmax(diff - tol))
            falhas.append(
                f"  coluna '{col}': {int(ruins.sum())} célula(s) fora da tolerância; "
                f"pior em linha {idx}: novo={a[idx]:.6e} baseline={b[idx]:.6e} "
                f"(|diff|={diff[idx]:.2e}, tol={tol[idx]:.2e})"
            )

    if falhas:
        print("Diferenças acima da tolerância (rtol=1e-5, atol=1e-8):")
        print("\n".join(falhas[:10]))
        return False

    return True


# ---------------------------------------------------------------------------
# Ponto de entrada
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    case_name = sys.argv[1]
    new_csv   = Path(sys.argv[2])

    # baseline sempre relativo ao diretório deste script
    script_dir   = Path(__file__).parent
    baseline_csv = script_dir / "baselines" / case_name / "RTSresult.csv"

    if not new_csv.exists():
        print(f"[ERRO] CSV não encontrado: {new_csv}")
        sys.exit(1)

    if not baseline_csv.exists():
        print(f"[ERRO] Baseline não encontrado: {baseline_csv}")
        sys.exit(1)

    ok = compare(case_name, new_csv, baseline_csv)

    if ok:
        print(f"[PASS] {case_name}")
        sys.exit(0)
    else:
        print(f"[FAIL] {case_name}")
        sys.exit(1)
