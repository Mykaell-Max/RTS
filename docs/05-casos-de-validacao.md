# 05 — Catálogo de Casos de Validação

> Inventário detalhado dos casos em [`validation/`](../validation), o que cada um testa,
> como rodar, e como vamos usá-los no trabalho de paralelização.
>
> Pré-requisitos: [00-comeca-aqui.md](00-comeca-aqui.md),
> [04-plano-de-ataque.md](04-plano-de-ataque.md).

---

## 1. Para Que Servem

Os casos em `validation/` são **benchmarks clássicos da literatura de transferência
radiativa**. Cada um reproduz um problema com resultado publicado em paper, e serve
como **oráculo numérico**: se nossa versão paralela gerar números diferentes do esperado,
sabemos que introduzimos bug.

São essenciais para o nosso trabalho por três razões:

1. **Regressão automática** — rodar após cada modificação para garantir que não quebrou
2. **Cobertura de caminhos** — cada caso exercita combinações diferentes de método/modelo
3. **Comparação com a tese** — os mesmos casos aparecem no Capítulo 5 da tese do Gustavo

---

## 2. Tabela Comparativa Rápida

| Caso | Dim | Malha padrão | Método | Modelo | Espalhamento | Gás | T não-uniforme |
|------|-----|-------------|--------|--------|--------------|-----|----------------|
| [`1D_Bordbar_WSGG`](#21-1d_bordbar_wsgg) | 1D | 30 | FAM | non-gray WSGG | — | CO₂/H₂O | sim |
| [`2D_Goutiere_flame`](#22-2d_goutiere_flame) | 2D | 21×21 | FAM | non-gray WSGG | — | CO₂/H₂O | sim |
| [`2D_Kim_scattering`](#23-2d_kim_scattering) | 2D | 25×25 | FAM | gray | anisotrópico linear | — | parede S=64.8K |
| [`2D_Shah_solution`](#24-2d_shah_solution) | 2D | 21×21 | FAM | gray | isotrópico | — | uniforme |
| [`3D_Bordbar_flame`](#25-3d_bordbar_flame) | 3D | 17×34×17 | FAM | non-gray WSGG | — | CO₂/H₂O | sim |
| **[`3D_Hsu_benchmark`](#26-3d_hsu_benchmark)** ⭐ | 3D | 21×21×21 | FAM | gray | isotrópico | — | uniforme, κ não-uniforme |
| [`3D_Soucasse_cavity`](#27-3d_soucasse_cavity) | 3D | 42×42×42 | FAM | gray | — | — | sim |
| [`demonstrations/symmetry_bc`](#28-demonstrationssymmetry_bc) | 3D | 20³ | FAM | gray | — | — | uniforme + simetria |

> ⭐ **`3D_Hsu_benchmark`** é o caso recomendado como benchmark principal de performance —
> 3D, FAM, absorção não-uniforme, tamanho ajustável.

---

## 3. Como Rodar um Caso (Workflow Atual)

O setup é **manual** hoje. Os passos típicos:

```bash
# 1. Copia o input específico do caso para a pasta input/
cp validation/3D_Hsu_benchmark/input.rts input/input.rts

# 2. Se o caso tem funções customizadas, substitui as do core
cp validation/3D_Hsu_benchmark/user_functions.f90 sources/RTS_functions.f90

# 3. Lê o readme — alguns pedem para editar output.rts manualmente
cat validation/3D_Hsu_benchmark/readme
# ... edita input/output.rts conforme instruído ...

# 4. Compila e roda
make clean && make run

# 5. Resultados aparecem em output/
ls output/
```

> ⚠️ **Esse processo é frágil:** alterações em `RTS_functions.f90` são esquecidas, o
> `output.rts` precisa ser editado de cabeça. Vamos automatizar isso na Fase 0 do plano.

### Script de automação proposto

Algo como `tests/run_validation.sh CASO`:

```bash
#!/bin/bash
# tests/run_validation.sh
CASE=$1
VALIDATION_DIR=validation/$CASE

# Backup do estado atual
cp input/input.rts input/input.rts.bak
cp sources/RTS_functions.f90 sources/RTS_functions.f90.bak

# Aplica config do caso
cp $VALIDATION_DIR/input.rts input/input.rts
if [ -f "$VALIDATION_DIR/user_functions.f90" ]; then
    cp $VALIDATION_DIR/user_functions.f90 sources/RTS_functions.f90
fi

# Compila e roda
make clean && make run

# Compara com referência (se existir)
if [ -d "tests/reference/$CASE" ]; then
    python tests/compare_fields.py output/ tests/reference/$CASE/
fi

# Restaura estado
mv input/input.rts.bak input/input.rts
mv sources/RTS_functions.f90.bak sources/RTS_functions.f90
```

---

## 4. Casos Detalhados

### 4.1 `1D_Bordbar_WSGG`

**Referência:** Bordbar et al. (2020), [DOI 10.1016/j.icheatmasstransfer.2019.104400](https://doi.org/10.1016/j.icheatmasstransfer.2019.104400)

**O que testa:** modelo **WSGG não-cinza** estendido para razões arbitrárias de CO₂/H₂O.

**Configuração:**
- **Dimensão:** 1D (mas com `nx=30, ny=40, nz=20` no input — só `x` é ativo)
- **Método:** FAM
- **Gás:** sim, com `XCO2=1e-4, XH2O=0.5`
- **Temperatura:** não-uniforme (via `user_functions.f90`)
- **Paredes:** West/East a 400K, demais a 300K
- **Absorção:** `κ = 10 m⁻¹`

**Por que importa para nós:**
- Único caso 1D — testa o caminho de código `DIMEN == 1`
- Exercita o `BAND_LOOP` (5 bandas WSGG)
- Validação numérica do `radiative_properties` (rotina WSGG)

**Custo computacional:** baixo (1D, malha pequena).

---

### 4.2 `2D_Goutiere_flame`

**Referência:** Goutiere et al. (2000), [DOI 10.1016/S0022-4073(99)00102-8](https://doi.org/10.1016/S0022-4073(99)00102-8)

**O que testa:** modelagem de gás real em recintos 2D.

**Configuração:**
- **Dimensão:** 2D (1m × 0.5m)
- **Método:** FAM
- **Gás:** sim (CO₂=0.1, H₂O=0.2)
- **Temperatura:** não-uniforme (via `user_functions.f90`)
- **Absorção:** `κ = 10 m⁻¹`
- **Paredes:** todas a 0K (idealização)

**Por que importa para nós:**
- Caminho 2D + WSGG (`DIMEN==2 && nongray_flag==.true.`)
- Testa as faces XY (não tem face Z não-trivial)

**Custo:** baixo (21×21, WSGG são 5 bandas).

---

### 4.3 `2D_Kim_scattering`

**Referência:** Kim & Lee (1988), [DOI 10.1016/0017-9310(88)90283-9](https://doi.org/10.1016/0017-9310(88)90283-9)

**O que testa:** **espalhamento anisotrópico** em recinto 2D.

**Configuração:**
- **Dimensão:** 2D (1m × 1m)
- **Método:** FAM
- **Gás:** não (caso acadêmico)
- **Absorção:** `κ = 0` (espalhamento puro)
- **Espalhamento:** `σ = 1 m⁻¹`, **anisotrópico**, fase linear F1
- **Parede South:** 64.8K, demais a 0K (irradia por uma face só)

**Por que importa para nós:**
- **Único caso com `aniso_flag = .true.`** — testa código de fase anisotrópica
- Testa o termo de espalhamento em `RHS_SM_FAM` (somatório duplo angular, O(N³·n²))
- Caminho mais caro do código

**Custo:** moderado (25² com espalhamento → loops angulares caros).

---

### 4.4 `2D_Shah_solution`

**Referência:** Shah, N.G. (1979) — Tese PhD, Imperial College London.
[Handle 10044/1/7839](http://hdl.handle.net/10044/1/7839)

**O que testa:** solução clássica de absorção/emissão em recinto 2D (sem espalhamento, sem gás).

**Configuração:**
- **Dimensão:** 2D (1m × 1m)
- **Método:** FAM
- **Gás:** não
- **Absorção:** `κ = 1 m⁻¹`
- **Temperatura uniforme:** 64.8K
- **Paredes:** todas a 0K

**Por que importa para nós:**
- Caso mais simples 2D — bom para sanity check
- Ideal para iniciar regressão (resultado conhecido analítico)

**Custo:** baixo.

---

### 4.5 `3D_Bordbar_flame`

**Referência:** Bordbar et al. (2014), [DOI 10.1016/j.combustflame.2014.03.013](https://doi.org/10.1016/j.combustflame.2014.03.013).
Caso originalmente proposto por Liu (1999).

**O que testa:** modelo WSGG 3D em geometria de chama (caixa alongada 2×4×2 m).

**Configuração:**
- **Dimensão:** 3D (`17×34×17`)
- **Método:** FAM
- **Gás:** sim (CO₂=0.85, H₂O=0.10 — oxi-combustão rica em CO₂)
- **Temperatura:** não-uniforme (via `user_functions.f90`)
- **Absorção:** `κ = 10 m⁻¹`
- **Paredes:** todas a 300K

**Por que importa para nós:**
- 3D + WSGG completo → caminho mais "produção" do código
- Exercita todos os 8 octantes
- Tamanho médio (~10k células) — viável para regressão razoavelmente rápida

**Custo:** médio (5 bandas × 8 octantes × 32 direções angulares).

---

### 4.6 `3D_Hsu_benchmark` ⭐

**Referência:** Hsu & Farmer (1997), [DOI 10.1115/1.2824087](https://doi.org/10.1115/1.2824087).

**O que testa:** cavidade 3D com **absorção não-uniforme** — comparação com solução Monte Carlo de referência.

**Configuração:**
- **Dimensão:** 3D cúbica (1m³)
- **Método:** FAM (mas DOM e P1 também funcionam — bom para testar todos)
- **Gás:** não
- **Temperatura uniforme:** 64.8K
- **Absorção e espalhamento:** **não-uniformes** (via `user_functions.f90`):
  ```fortran
  κ(x,y,z) = 0.9·(1-2|x-0.5|)·(1-2|y-0.5|)·(1-2|z-0.5|) + 0.1
  σ(x,y,z) = 0.9·κ(x,y,z)
  ```
  Resultado: bolha de absorção/espalhamento no centro, decaindo para as paredes.
- **Paredes:** todas a 0K

**Por que importa para nós (e por que recomendamos como benchmark principal):**
- 3D real, exercitando todos os 8 octantes
- Absorção e espalhamento **espacialmente variáveis** → impossível otimizar com simplificações
- Resultado Monte Carlo conhecido com **alta precisão** (dá oráculo confiável)
- Geometria simples — fácil escalar o tamanho da malha (21³ → 64³ → 128³ → 256³) mantendo o problema físico igual

**Custo:** ajustável conforme malha:

| nx=ny=nz | Tempo serial estimado | Memória `IG` (FAM 32 dir) |
|----------|----------------------|---------------------------|
| 21 (default) | ~segundos | ~3 MB |
| 64 | ~minutos | ~70 MB |
| 128 | dezenas de min | ~550 MB |
| 256 | horas | ~4 GB |

---

### 4.7 `3D_Soucasse_cavity`

**Referência:** Soucasse et al. (2014), [DOI 10.1615/ComputThermalScien.2012005118](https://doi.org/10.1615/ComputThermalScien.2012005118)

**O que testa:** acoplamento de radiação molecular com convecção natural em cavidade cúbica aquecida diferencialmente.

**Configuração:**
- **Dimensão:** 3D cúbica (1m³)
- **Malha:** 42×42×42 (já vem maior — mais demorada)
- **Método:** FAM
- **Gás:** não (gray, com `κ = 1`)
- **Temperatura:** não-uniforme (via `user_functions.f90`)
- **Paredes:** todas a 300K com **emissividade 0.5** (paredes não-negras!)

**Por que importa para nós:**
- Caso com **paredes não-negras** (emissividade ≠ 1) — testa `epsilon_rad` no código
- Malha já razoável (42³ ≈ 74k células) — bom benchmark intermediário sem precisar escalar

**Custo:** médio-alto.

---

### 4.8 `demonstrations/symmetry_bc`

**Referência:** demonstração interna (sem paper)

**O que testa:** **condições de contorno de simetria** em paredes.

**Configuração:**
- **Dimensão:** 3D (0.5m³, malha 20³)
- **Método:** FAM
- **Gás:** não
- **Temperatura uniforme:** 64.8K
- **Paredes:** East, South, Bottom com **`SBCwall = .true.`** (simetria)

**Por que importa para nós:**
- Único caso com BCs de simetria — testa o caminho do código para `SBCwall`
- Sanity check de comportamento físico (simetria = espelho)

**Custo:** baixo.

---

## 5. Cobertura de Caminhos do Código

Mapeamento "qual caso exercita qual funcionalidade":

| Funcionalidade | Casos que exercitam |
|---------------|---------------------|
| `DIMEN == 1` | `1D_Bordbar_WSGG` |
| `DIMEN == 2` | `2D_Kim_scattering`, `2D_Shah_solution`, `2D_Goutiere_flame` |
| `DIMEN == 3` | `3D_Hsu_benchmark`, `3D_Bordbar_flame`, `3D_Soucasse_cavity`, `symmetry_bc` |
| `nongray_flag = .true.` (WSGG) | `1D_Bordbar_WSGG`, `2D_Goutiere_flame`, `3D_Bordbar_flame` |
| `nongray_flag = .false.` (gray) | `2D_Kim`, `2D_Shah`, `3D_Hsu`, `3D_Soucasse`, `symmetry_bc` |
| `aniso_flag = .true.` (espalhamento) | `2D_Kim_scattering` |
| `gas_prop = .true.` | todos os WSGG (`1D_Bordbar`, `2D_Goutiere`, `3D_Bordbar`) |
| BCs simetria | `symmetry_bc` |
| Paredes não-negras | `3D_Soucasse_cavity` (ε=0.5) |
| Campo `κ` não-uniforme | `3D_Hsu_benchmark` |
| Campo T não-uniforme | quase todos com `user_functions.f90` |

> **Conclusão:** rodar todos os 8 casos cobre virtualmente todos os caminhos de código.
> Para regressão rápida durante desenvolvimento, podemos rodar um subset
> (ex.: `2D_Shah` + `3D_Hsu`). Para validação completa antes de aceitar fase, rodar tudo.

---

## 6. Como Vamos Usar Cada Caso por Fase

| Fase do plano | Casos usados | Propósito |
|---------------|--------------|-----------|
| **F0 baseline** | todos | Medir tempos, gerar arquivos de referência |
| **F1 OpenMP fácil** | `3D_Hsu_benchmark` (rápido), depois todos antes de fechar fase | Regressão rápida no dev, completa antes de merge |
| **F2 OpenMP difícil** | todos (especialmente `3D_Bordbar_flame` para SOR via energia) | SOR e wavefront podem ter pequenas diferenças numéricas — comparar tolerância |
| **F3 refactor local** | todos | Mudança estrutural — testar tudo |
| **F4 MPI** | `3D_Hsu_benchmark` em vários tamanhos + todos antes de fechar | Validar com 1, 2, 4, 8, 16 ranks |
| **F5 híbrido** | `3D_Hsu_benchmark` grande (128³+) | Benchmark de escala |
| **F6 benchmarks** | todos + variações de tamanho | Tabela final |

---

## 7. Estratégia de Comparação Numérica

Para regressão automática precisamos definir **o que é "resultado igual"**:

### 7.1 Campos a comparar

Cada caso produz na saída (em `output/`):
- `G(i,j,k)` — incidência radiativa
- `S_rad(i,j,k)` — fonte radiativa
- `Q_radw(i,j,k)` — fluxos nas paredes

### 7.2 Métricas

Para cada campo, comparar com referência usando:

- **Erro máximo absoluto:** `max(|atual - ref|)`
- **Erro RMS:** `sqrt(mean((atual - ref)²))`
- **Erro relativo:** `max(|atual - ref| / |ref|)` onde `|ref| > tolerância`

### 7.3 Tolerâncias

| Tipo de mudança | Tolerância sugerida |
|----------------|---------------------|
| OpenMP em loops sem dependência | **Bit-a-bit idêntico** (mesma ordem de operações) |
| OpenMP em `BAND_LOOP` com REDUCTION | Erro relativo ≤ 1e-12 (ordem de soma muda) |
| Red-black SOR | Erro relativo ≤ 1e-6 (convergência diferente) |
| MPI sem mudança algorítmica | Erro relativo ≤ 1e-10 |
| MPI com KBA wavefront | Erro relativo ≤ 1e-8 |

Acima da tolerância → **falha de regressão**, investigar.

---

## 8. Action Items Imediatos (Fase 0)

Para integrar essa estratégia ao trabalho:

1. **Criar `tests/run_validation.sh`** que automatiza o setup de caso
2. **Criar `tests/compare_fields.py`** que compara dois conjuntos de output `.dat` com tolerância configurável
3. **Rodar todos os 8 casos em sequencial** e gerar `tests/reference/CASO/` com os outputs de referência
4. **Documentar tempos serial** de cada caso no `docs/05-baseline.md` (futuro)
5. **Adicionar `make validate`** ao Makefile que roda todos os casos sequencialmente e reporta resultado

---

## 9. Limitações dos Casos para Benchmark de Performance

Apesar de excelentes para validação, todos os casos vêm com **malhas pequenas** (21³ típico).
Isso é ótimo para regressão rápida, mas **insuficiente para benchmark de paralelização** —
um caso que roda em 2 segundos serial não mostra speedup claro com 16 cores.

**Estratégia:**

- **Validação numérica:** rodar com malhas padrão (rápido, conclusivo)
- **Benchmark de performance:** rodar `3D_Hsu_benchmark` em tamanhos crescentes
  (21³ → 64³ → 128³ → 256³), fixando o resto da física

> O caso `3D_Hsu_benchmark` é único nesse aspecto: como a geometria é uma cavidade
> 1m³ cúbica e os campos são definidos por funções contínuas (`user_functions.f90`),
> escalar a malha é trivial — só mudar `nx,ny,nz` no `input.rts`. O problema físico
> permanece o mesmo, só com mais resolução.

---

*Documento será atualizado conforme novos casos forem adicionados ou tolerâncias forem
calibradas com base em experiência.*
