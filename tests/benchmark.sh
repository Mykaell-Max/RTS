#!/usr/bin/env bash
# =============================================================================
# benchmark.sh  --  valida a Fase 2 e mede speedup no caso mais pesado
#
# Uso (da raiz do repositório):
#     bash tests/benchmark.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"

# cores
if [ -t 1 ]; then
    GREEN="\033[0;32m"; RED="\033[0;31m"; YELLOW="\033[0;33m"
    RESET="\033[0m";    BOLD="\033[1m"
else
    GREEN=""; RED=""; YELLOW=""; RESET=""; BOLD=""
fi

echo ""
echo -e "${BOLD}=============================================${RESET}"
echo -e "${BOLD}  RTS  --  Regressão + Benchmark Fase 2     ${RESET}"
echo -e "${BOLD}=============================================${RESET}"
echo ""

# --- 1. Regressão com build debug (fcheck=all) ---
echo -e "${BOLD}>>> Passo 1: regressão OpenMP (build debug, 8 threads)${RESET}"
export OMP_NUM_THREADS=8
bash tests/run_all.sh
echo ""

# --- 2. Build otimizado para benchmark ---
echo -e "${BOLD}>>> Passo 2: build otimizado (-O2) para benchmark de speedup${RESET}"
CASE_RTS="$REPO_ROOT/validation/3D_Soucasse/RTS"
cp "$CASE_RTS/input/"*.rts input/
cp "$CASE_RTS/sources/RTS_functions.f90" sources/
make clean -s
make OMP="-fopenmp" FFLAGS="-c -O2 -fopenmp" -s
echo "Build concluido."
echo ""

# --- 3. Speedup ---
echo -e "${BOLD}>>> Passo 3: speedup por numero de threads (3D_Soucasse, -O2)${RESET}"
echo ""
printf "%-12s  %s\n" "Threads" "Tempo real"
printf "%-12s  %s\n" "-------" "----------"

T1=""
for N in 1 2 4 8 16; do
    export OMP_NUM_THREADS=$N
    TEMPO=$( { time ./bin/RTS > /dev/null; } 2>&1 | grep real | awk '{print $2}' )
    SEGS=$(echo "$TEMPO" | sed 's/m/:/g' | awk -F: '{printf "%.3f", $1*60 + $2}')
    if [ -z "$T1" ]; then T1="$SEGS"; fi
    if [ "$T1" != "0" ] && [ -n "$T1" ]; then
        SPEEDUP=$(awk "BEGIN {printf \"%.2f\", $T1/$SEGS}")
    else
        SPEEDUP="1.00"
    fi
    printf "%-12s  %s  (speedup: %sx)\n" "$N" "$TEMPO" "$SPEEDUP"
done

echo ""
echo -e "${BOLD}=============================================${RESET}"
echo -e "${GREEN}Benchmark concluido.${RESET}"
echo ""
