#!/usr/bin/env bash
# =============================================================================
# run_all.sh  --  compila e valida todos os 7 casos de regressão do RTS
#
# Uso (da raiz do repositório):
#     bash tests/run_all.sh                  # compara contra baselines salvas
#     bash tests/run_all.sh --update-baselines  # só na 1ª vez: gera as baselines
#
# Fluxo por caso:
#   1. Copia input/*.rts e RTS_functions.f90 do caso para o root
#   2. make  (recompila apenas o que mudou)
#   3. ./bin/RTS  (gera output/RTSresult.csv)
#   4. check_regression.py  (compara com tests/baselines/CASO/RTSresult.csv)
#      --ou--  copia resultado para tests/baselines/ se --update-baselines
#
# Pré-requisito: compilação inicial com  make clean && make  antes da 1ª rodada
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

UPDATE_MODE=0
if [[ "${1:-}" == "--update-baselines" ]]; then
    UPDATE_MODE=1
fi

CASOS=(
    1D_Bordbar
    2D_Goutiere
    2D_Kim
    2D_Shah
    3D_Bordbar
    3D_Hsu
    3D_Soucasse
)

# cores (desativa se não for terminal)
if [ -t 1 ]; then
    GREEN="\033[0;32m"; RED="\033[0;31m"; YELLOW="\033[0;33m"
    RESET="\033[0m";    BOLD="\033[1m"
else
    GREEN=""; RED=""; YELLOW=""; RESET=""; BOLD=""
fi

declare -A RESULTADOS
FALHAS=0

# Modo de build:
#   - baseline (--update-baselines): binário SERIAL (OMP=) -> referência confiável
#   - regressão (padrão):            binário OpenMP        -> usa OMP_NUM_THREADS
if [ "$UPDATE_MODE" -eq 1 ]; then
    MAKE_OMP="OMP="          # serial
    BUILD_DESC="serial (referência)"
else
    MAKE_OMP=""              # default do makefile = -fopenmp
    BUILD_DESC="OpenMP (OMP_NUM_THREADS=${OMP_NUM_THREADS:-padrão do SO})"
fi

echo ""
if [ "$UPDATE_MODE" -eq 1 ]; then
    echo -e "${BOLD}=============================================${RESET}"
    echo -e "${BOLD}  RTS  --  Gerando baselines (7 casos)       ${RESET}"
    echo -e "${BOLD}=============================================${RESET}"
else
    echo -e "${BOLD}=============================================${RESET}"
    echo -e "${BOLD}  RTS  --  Regressão completa (7 casos)      ${RESET}"
    echo -e "${BOLD}=============================================${RESET}"
fi
echo -e "  build: ${BUILD_DESC}"
echo ""

# Limpeza única no início: garante que os .o correspondem ao modo de build
# atual (serial vs OpenMP). Dentro de uma rodada o flag é constante, então o
# make incremental por caso (que só recompila RTS_functions) é suficiente.
(cd "$REPO_ROOT" && make clean -s $MAKE_OMP) 2>/dev/null || true

for CASO in "${CASOS[@]}"; do
    CASE_RTS="$REPO_ROOT/validation/$CASO/RTS"
    NEW_CSV="$REPO_ROOT/output/RTSresult.csv"
    BASELINE="$REPO_ROOT/tests/baselines/$CASO/RTSresult.csv"

    echo -e "${BOLD}--- $CASO ---${RESET}"

    # --- 1. Copiar inputs e RTS_functions do caso ---
    cp "$CASE_RTS/input/"*.rts "$REPO_ROOT/input/"
    cp "$CASE_RTS/sources/RTS_functions.f90" "$REPO_ROOT/sources/RTS_functions.f90"

    # --- 2. Compilar (incremental: só recompila o que mudou) ---
    if ! (cd "$REPO_ROOT" && make -s $MAKE_OMP) ; then
        echo -e "  ${RED}[FAIL] compilação falhou${RESET}"
        RESULTADOS[$CASO]="FAIL (compilação)"
        FALHAS=$((FALHAS + 1))
        continue
    fi

    # --- 3. Executar ---
    (cd "$REPO_ROOT" && make cleanout -s)
    if ! (cd "$REPO_ROOT" && ./bin/RTS > /dev/null 2>&1) ; then
        echo -e "  ${RED}[FAIL] execução falhou${RESET}"
        RESULTADOS[$CASO]="FAIL (execução)"
        FALHAS=$((FALHAS + 1))
        continue
    fi

    # --- 4. Comparar ou atualizar baseline ---
    if [ "$UPDATE_MODE" -eq 1 ]; then
        mkdir -p "$REPO_ROOT/tests/baselines/$CASO"
        cp "$NEW_CSV" "$BASELINE"
        echo -e "  ${YELLOW}[UPDATED] baseline salva${RESET}"
        RESULTADOS[$CASO]="UPDATED"
    else
        if python3 "$SCRIPT_DIR/check_regression.py" "$CASO" "$NEW_CSV" ; then
            RESULTADOS[$CASO]="PASS"
        else
            RESULTADOS[$CASO]="FAIL (regressão)"
            FALHAS=$((FALHAS + 1))
        fi
    fi

    echo ""
done

# --- sumário final ---
echo ""
echo -e "${BOLD}=============================================${RESET}"
echo -e "${BOLD}  Resultado                                  ${RESET}"
echo -e "${BOLD}=============================================${RESET}"
for CASO in "${CASOS[@]}"; do
    RES="${RESULTADOS[$CASO]}"
    if [ "$RES" = "PASS" ]; then
        echo -e "  ${GREEN}PASS   ${RESET} $CASO"
    elif [ "$RES" = "UPDATED" ]; then
        echo -e "  ${YELLOW}UPDATED${RESET} $CASO"
    else
        echo -e "  ${RED}FAIL   ${RESET} $CASO  -- $RES"
    fi
done
echo ""

if [ "$UPDATE_MODE" -eq 1 ]; then
    echo -e "${YELLOW}Baselines atualizadas. Rode sem --update-baselines para validar próximas mudanças.${RESET}"
    exit 0
elif [ "$FALHAS" -eq 0 ]; then
    echo -e "${GREEN}Todos os casos aprovados.${RESET}"
    exit 0
else
    echo -e "${RED}$FALHAS caso(s) reprovado(s).${RESET}"
    exit 1
fi
