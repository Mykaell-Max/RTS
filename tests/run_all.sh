#!/usr/bin/env bash
# =============================================================================
# run_all.sh  --  compila e valida todos os 7 casos de regressão do RTS
#
# Uso (da raiz do repositório):
#     bash tests/run_all.sh
#
# O que faz por caso:
#   1. cd validation/CASO/RTS
#   2. make clean && make        (compila com as sources embarcadas do caso)
#   3. ./bin/RTS                 (executa, gera output/RTSresult.csv)
#   4. check_regression.py       (compara com tests/baselines/CASO/RTSresult.csv)
#
# Resultado final: tabela PASS/FAIL + exit 0 se todos passaram, 1 se algum falhou
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

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
    GREEN="\033[0;32m"; RED="\033[0;31m"; RESET="\033[0m"; BOLD="\033[1m"
else
    GREEN=""; RED=""; RESET=""; BOLD=""
fi

declare -A RESULTADOS
FALHAS=0

echo ""
echo -e "${BOLD}=============================================${RESET}"
echo -e "${BOLD}  RTS  --  Regressão completa (7 casos)     ${RESET}"
echo -e "${BOLD}=============================================${RESET}"
echo ""

for CASO in "${CASOS[@]}"; do
    CASE_DIR="$REPO_ROOT/validation/$CASO/RTS"
    NEW_CSV="$CASE_DIR/output/RTSresult.csv"

    echo -e "${BOLD}--- $CASO ---${RESET}"

    # 1. Compilar
    if ! (cd "$CASE_DIR" && make clean -s && make -s) ; then
        echo -e "  ${RED}[FAIL] compilação falhou${RESET}"
        RESULTADOS[$CASO]="FAIL (compilação)"
        FALHAS=$((FALHAS + 1))
        continue
    fi

    # 2. Executar
    if ! (cd "$CASE_DIR" && ./bin/RTS > /dev/null 2>&1) ; then
        echo -e "  ${RED}[FAIL] execução falhou${RESET}"
        RESULTADOS[$CASO]="FAIL (execução)"
        FALHAS=$((FALHAS + 1))
        continue
    fi

    # 3. Comparar contra baseline
    if python3 "$SCRIPT_DIR/check_regression.py" "$CASO" "$NEW_CSV" ; then
        RESULTADOS[$CASO]="PASS"
    else
        RESULTADOS[$CASO]="FAIL (regressão)"
        FALHAS=$((FALHAS + 1))
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
        echo -e "  ${GREEN}PASS${RESET}  $CASO"
    else
        echo -e "  ${RED}FAIL${RESET}  $CASO  ($RES)"
    fi
done
echo ""

if [ "$FALHAS" -eq 0 ]; then
    echo -e "${GREEN}Todos os casos aprovados.${RESET}"
    exit 0
else
    echo -e "${RED}$FALHAS caso(s) reprovado(s).${RESET}"
    exit 1
fi
