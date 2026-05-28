# 04 — Plano de Ataque (Paralelização do RTS standalone)

> **Documento de execução.** Aqui está o roadmap concreto, dividido em fases, para
> paralelizar o solver do RTS. Cada fase tem **objetivo, arquivos afetados, critério de
> aceitação, esforço estimado e riscos**.
>
> Pré-requisitos de leitura: [00-comeca-aqui.md](00-comeca-aqui.md),
> [01-arquitetura-atual.md](01-arquitetura-atual.md),
> [02-relatorio-mfsim-mpi.md](02-relatorio-mfsim-mpi.md).

---

## 1. Escopo

### 1.1 O que está incluso

- Refatoração dos arquivos `sources/RTS_*.f90` para suportar paralelismo
- Adição de OpenMP nos loops apropriados
- Adição de MPI com decomposição espacial (KBA/wavefront)
- Validação numérica contra os casos em `validation/`
- Instrumentação e medição de performance

### 1.2 O que NÃO está incluso (fora de escopo)

| Item | Por quê |
|------|---------|
| Refatorar `rad_interface.f90` ("blocão") | Repositório MFSim — papel da equipe deles |
| Refatorar `rad_core.f90`, `RTS_connection.f90` | Repositório MFSim — papel da equipe deles |
| Implementar RadCal, FSCK, SNB | Fora do escopo (modelos de absorção) |
| Versão GPU (CUDA/OpenACC) | Fase muito posterior, se aplicável |
| Suporte a malhas não-estruturadas | Fora do escopo |

### 1.3 Pressupostos

- Temos acesso de escrita ao código em `sources/`
- Temos pelo menos uma máquina multi-core para desenvolvimento
- Para validar MPI: precisamos eventualmente de cluster ou máquina com vários sockets
- Não dependemos de avanços no MFSim para começar — fazemos o RTS "MPI-ready"

---

## 2. Os 3 Problemas Que Vamos Atacar

Recap dos problemas que sobram após excluir o "blocão" (Frente B):

| Problema | Causa estrutural | Como atacamos |
|----------|------------------|---------------|
| 🐢 **Tempo do solver** | Loops `O(N³ · n_ang²)` em série | OpenMP nos loops independentes + MPI nos sweeps |
| 💾 **Memória** | `IG(nxi,nyi,nzi,nt,np)` alocado por inteiro em cada rank | Decomposição espacial: cada rank aloca só sua fatia |
| 📉 **Não escala** | Solver é todo serial | OpenMP intra-nó + MPI inter-nó |

---

## 3. Estratégia em Camadas

```mermaid
flowchart TB
    subgraph L1["🟢 Camada 1 — OpenMP (intra-nó)"]
        direction TB
        L1A[Loops embaraçosamente paralelos]
        L1B[Loops com dependência: red-black, wavefront]
    end

    subgraph L2["🔵 Camada 2 — MPI (inter-nó)"]
        direction TB
        L2A[Decomposição espacial cartesiana]
        L2B[Halo exchange]
        L2C[KBA nos sweeps]
        L2D[Reduções globais para resíduo]
    end

    subgraph L3["🟣 Camada 3 — Híbrido]"]
        L3A[MPI inter-nó + OpenMP intra-nó]
    end

    L1 --> L2
    L2 --> L3
```

**Por que essa ordem:**

1. **OpenMP primeiro** porque traz ganho rápido mesmo sem mexer na estrutura. É a "rede de segurança" — se o MPI atrasar, ainda temos speedup.
2. **MPI depois** porque exige refatoração estrutural (arrays locais, halos, mudança nos sweeps). Mas é o que resolve o problema de memória e escala real.
3. **Híbrido no final** porque só faz sentido depois das duas camadas funcionarem isoladas.

---

## 4. Plano por Fases

> **Estimativas de esforço** são em "semanas-equivalentes" assumindo 1 pessoa em dedicação
> moderada. Multiplicar conforme o time real.

### Fase 0 — Baseline e Instrumentação (foundation)

**Objetivo:** estabelecer pontos de medição antes de mexer em qualquer coisa.
**Esforço:** ~1 semana
**Pré-requisito:** ambiente de build (gfortran + make) funcionando

#### Atividades

1. **Compilar e rodar** o RTS atual (serial) com um caso de validação
2. **Escolher caso de benchmark** principal (sugestão: `3D_Hsu_benchmark` — 3D, gas absorvente, suficientemente caro para diferenças aparecerem)
3. **Adicionar timers** nas rotinas principais:
   ```fortran
   double precision :: t0, t1
   call cpu_time(t0)
   ! ... rotina ...
   call cpu_time(t1)
   write(*,*) 'tempo XYZ:', t1-t0
   ```
   Pontos a instrumentar: `radiative_properties`, `BAND_LOOP`, `RHS_SM_DOM/FAM`,
   `orthogonal_loop`/`agular_loop`, `SOR`, total do `rtesolve`
4. **Salvar resultados de referência** (campos `G`, `S_rad`, `Q_radw`) em arquivos `.dat`
   binários para comparação byte-a-byte depois
5. **Criar script de regressão** — roda o caso, compara saída atual com referência, falha se erro > tolerância numérica
6. **Documentar baseline** — tabela: tempo total, tempo por rotina, em N=64, N=128, N=256

#### Arquivos afetados

| Arquivo | Mudança |
|---------|---------|
| [`RTS_radiation.f90`](sources/RTS_radiation.f90) | Adicionar `cpu_time` em pontos-chave |
| [`RTS_solvers.f90`](sources/RTS_solvers.f90) | Adicionar `cpu_time` no `SOR` |
| Novo: `tests/regression.sh` | Script de comparação |
| Novo: `docs/baseline.md` | Tabela de medições |

#### Critério de aceitação

- [ ] RTS compila e roda o caso de validação escolhido
- [ ] Tabela com tempo de cada rotina principal
- [ ] Script de regressão funcionando (detecta se mudança quebra o numérico)
- [ ] Arquivos de referência salvos em `tests/reference/`

#### Riscos

- Caso 3D escolhido pode ser pequeno demais (rodar segundos) → escolher tamanho que rode em ~minutos para ter sinal claro
- Diferenças entre máquinas (gfortran 9 vs 10, flags) → fixar versão e flags no Makefile

---

### Fase 1 — OpenMP nos Loops Embaraçosamente Paralelos

**Objetivo:** primeiro speedup, sem mexer em estrutura.
**Esforço:** ~2 semanas
**Pré-requisito:** Fase 0 completa

#### Atividades

1. **Habilitar OpenMP no Makefile**
   ```makefile
   FFLAGS = ... -fopenmp
   # No linker também:
   $(FC) -fopenmp ...
   ```

2. **Paralelizar `RHS_SM_DOM` e `RHS_SM_FAM`** — loops 5D sem dependência:
   ```fortran
   !$OMP PARALLEL DO PRIVATE(i,j,l,ls,m,SMSUM,aux) SCHEDULE(STATIC)
   do k=2,nzp
       ...
   end do
   !$OMP END PARALLEL DO
   ```

3. **Paralelizar `radiative_properties` (WSGG)** — triplo loop espacial célula a célula

4. **Paralelizar `G_DOM`, `G_FAM`** — soma angular por célula, independente entre células

5. **Paralelizar `wall_fluxes_DOM`, `wall_fluxes_FAM`, `wall_fluxes_pone`** — loops 2D nas faces

6. **Paralelizar `BAND_LOOP`** (mais delicado por causa de redução):
   ```fortran
   !$OMP PARALLEL DO PRIVATE(IBND,IGband,Gband,Sband,IBFF) &
   !$OMP&            REDUCTION(+:G,S_rad) SCHEDULE(STATIC)
   BAND_LOOP: do IBND = 1,nsb
       ...
   end do BAND_LOOP
   !$OMP END PARALLEL DO
   ```
   Atenção: `cappa`, `beta` precisam ser tornadas privadas ou refatoradas para arrays
   indexados por banda

7. **Rodar regressão** após cada modificação

#### Arquivos afetados

| Arquivo | Mudança |
|---------|---------|
| [`makefile`](makefile) | Adicionar `-fopenmp` |
| [`RTS_radiation.f90`](sources/RTS_radiation.f90) | Diretivas OpenMP nas rotinas listadas |
| [`RTS_absorption.f90`](sources/RTS_absorption.f90) | Diretivas OpenMP em `radiative_properties` |

#### Critério de aceitação

- [ ] Compila com `-fopenmp` sem warnings de race condition
- [ ] Speedup ≥ **4×** em 8 threads (em FAM 3D non-gray, caso de benchmark)
- [ ] Speedup ≥ **8×** em 16 threads
- [ ] Resultado numérico **bit-a-bit idêntico** ao serial (a maioria das paralelizações
      aqui é determinística; se diferir, pode ser ordem de soma — verificar)
- [ ] Script de regressão verde

#### Riscos

| Risco | Mitigação |
|-------|-----------|
| Race condition em `cappa`/`beta` no `BAND_LOOP` | Tornar locais à iteração ou usar `firstprivate` |
| Speedup pobre por false sharing | Verificar layout de memória, usar `SCHEDULE(STATIC, chunk)` adequado |
| Resultado numérico difere por ordem de soma | Aceitável se for diferença na 14ª casa; se for maior, investigar |

---

### Fase 2 — OpenMP nos Loops Difíceis (sweeps e SOR)

**Objetivo:** atacar os hot loops que têm dependência interna.
**Esforço:** ~3-4 semanas
**Pré-requisito:** Fase 1 funcionando

#### Atividades

1. **Red-black ordering no SOR**
   - Trocar Gauss-Seidel sequencial por padrão "vermelho-preto" (xadrez 3D)
   - Células "vermelhas" (i+j+k par) atualizam usando "pretas"; depois invertem
   - Cada fase é embaraçosamente paralela
   - **Implica mudança numérica** (convergência diferente, pode precisar de mais iterações)

2. **Wavefront diagonal nos sweeps DOM/FAM**
   - Em vez de varrer `(k,j,i)` sequencial, varrer **diagonais** `i+j+k = const`
   - Todas as células da mesma diagonal são independentes
   - Implementação típica: loop externo na diagonal, loop interno paralelizado
   ```fortran
   ! pseudocódigo octante 1
   do diag = 0, (nxp-2)+(nyp-2)+(nzp-2)
     !$OMP PARALLEL DO PRIVATE(i,j,k,...)
     do k = max(2,diag-nxp-nyp+4), min(nzp, diag+2)
        do j = ...
           i = diag - j - k + 6
           if (i in [2,nxp]) then
              ! update IG(i,j,k,l,m)
           end if
        end do
     end do
     !$OMP END PARALLEL DO
   end do
   ```
   - Repetir para cada octante (8×) e cada direção angular (l,m)

3. **Alternativa simples no sweep:** paralelizar só a direção angular `(l,m)` por célula
   - Mais simples mas potencial menor — fica como plano B

#### Arquivos afetados

| Arquivo | Mudança |
|---------|---------|
| [`RTS_solvers.f90`](sources/RTS_solvers.f90) | Refatorar `SOR` para red-black |
| [`RTS_radiation.f90`](sources/RTS_radiation.f90) | Refatorar `orthogonal_loop` e `agular_loop` para wavefront |

#### Critério de aceitação

- [ ] Convergência do SOR similar à atual (talvez 1.2× mais iterações é aceitável)
- [ ] Speedup adicional ≥ **2×** em cima da Fase 1
- [ ] Validação numérica passa em todos os casos de `validation/`

#### Riscos

| Risco | Mitigação |
|-------|-----------|
| Red-black piora convergência demais | Ter SOR sequencial como fallback (flag em runtime) |
| Wavefront é complexo de implementar certo | Começar com versão 2D, testar, depois 3D |
| Validação física diverge | Comparar com benchmark serial: erro absoluto deve estar dentro de `rad_tol` |

---

### Fase 3 — Preparação Estrutural para MPI

**Objetivo:** refatorar o código para que cada rank possa eventualmente trabalhar só com
sua fatia. **Ainda não usa MPI** — só prepara o terreno.
**Esforço:** ~3 semanas
**Pré-requisito:** Fase 1 (Fase 2 opcional)

#### Atividades

1. **Introduzir conceito de "subdomínio local" no `RTS_global.f90`**
   ```fortran
   ! adicionar:
   integer :: nxi_loc, nyi_loc, nzi_loc       ! tamanhos locais
   integer :: i_offset, j_offset, k_offset    ! deslocamento global do meu pedaço
   integer :: my_rank, n_ranks                ! identidade MPI
   integer :: nbr_west, nbr_east, nbr_south, nbr_north, nbr_bottom, nbr_top
   ```
   No serial atual, `nxi_loc = nxi`, `offset = 0`, `n_ranks = 1`, vizinhos = `MPI_PROC_NULL`.
   Isso permite o código compilar e rodar **igual ao atual** mesmo após a refatoração.

2. **Substituir todas as alocações** `(nxi, nyi, nzi)` por `(nxi_loc, nyi_loc, nzi_loc)`
   nas centenas de declarações

3. **Distinguir BC física vs interna no `RTS_bc.f90`**
   - Adicionar flags `is_physical_west`, `is_physical_east`, etc.
   - Hoje todas são `.true.` (caso serial)
   - Quando virar MPI, faces internas serão `.false.` e usarão halo

4. **Refatorar `ETANORM`, `MAXVAL`** para terem versão "global-ready"
   - Pode ser uma macro/função que no serial é `MAXVAL(...)` e no MPI vira `MPI_ALLREDUCE`

5. **Mover I/O para "root-only"**
   - `input_data` em `RTS_input.f90`: só rank 0 lê o arquivo, depois `MPI_Bcast` (no MPI)
   - No serial isso fica como está

#### Arquivos afetados

Quase todos:
- [`RTS_global.f90`](sources/RTS_global.f90)
- [`RTS_start.f90`](sources/RTS_start.f90)
- [`RTS_radiation.f90`](sources/RTS_radiation.f90)
- [`RTS_solvers.f90`](sources/RTS_solvers.f90)
- [`RTS_bc.f90`](sources/RTS_bc.f90)
- [`RTS_energy.f90`](sources/RTS_energy.f90)

#### Critério de aceitação

- [ ] Código compila e roda com mesma performance da Fase 2
- [ ] Resultado numérico idêntico (regressão verde)
- [ ] Variáveis de subdomínio existem (mesmo que com valores triviais)
- [ ] BC física vs interna está abstrata, mesmo que sempre física no serial

#### Riscos

| Risco | Mitigação |
|-------|-----------|
| Refatoração massiva quebra coisas sutis | Commits pequenos, rodar regressão após cada |
| Nomes de variáveis ficam confusos | Convenção clara: sufixo `_loc` para local |

---

### Fase 4 — MPI com Decomposição Espacial (o trabalho pesado)

**Objetivo:** RTS roda em paralelo de verdade, escala em cluster.
**Esforço:** ~6-8 semanas
**Pré-requisito:** Fase 3 obrigatória

#### Atividades

1. **`MPI_Init`/`MPI_Finalize` em `RTS_main.f90`**
   - Em modo standalone, o RTS inicializa MPI
   - Em modo acoplado ao MFSim, MFSim já inicializou — usar `MPI_Initialized` para detectar

2. **Topologia cartesiana** em `RTS_start.f90`
   ```fortran
   call MPI_Cart_create(MPI_COMM_WORLD, 3, [Px,Py,Pz], &
                        [.false.,.false.,.false.], .true., RTS_COMM, ierr)
   call MPI_Cart_coords(RTS_COMM, my_rank, 3, my_coords, ierr)
   call MPI_Cart_shift(RTS_COMM, 0, 1, nbr_west, nbr_east, ierr)
   ! ... etc para y e z
   ```

3. **Calcular subdomínio de cada rank**
   ```fortran
   nx_loc = nx / Px + (1 if my_coords(1) < mod(nx, Px) else 0)
   i_offset = ... ! deslocamento global
   nxi_loc = nx_loc + 2  ! +2 para halo
   ```

4. **Halo exchange para campos auxiliares** (`T_energy`, `cappa`, `sigma`, ...)
   - Após inicialização e a cada update
   - 6 trocas (oeste, leste, sul, norte, fundo, topo)
   - `MPI_Sendrecv` com buffer das faces

5. **KBA nos sweeps DOM/FAM** — o core do trabalho
   - Para cada octante, identificar quais 3 faces são "entrada" (upstream) e quais 3 são "saída"
   - `MPI_Irecv` posts antes do sweep
   - `MPI_Wait` antes de processar células próximas à face
   - `MPI_Isend` após processar células próximas à face de saída
   - Iterar sobre as 8 octantes (cada um define ordem de varredura diferente)
   - Pseudocódigo detalhado em [02-relatorio-mfsim-mpi.md §5.3](02-relatorio-mfsim-mpi.md)

6. **Resíduo global**
   ```fortran
   call MPI_ALLREDUCE(eps_local, eps_global, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, RTS_COMM, ierr)
   if (eps_global < rad_tol) exit
   ```

7. **SOR distribuído**
   - Cada iteração: halo exchange + computação local + ALLREDUCE no resíduo
   - Pode ficar lento (muitas comunicações pequenas) — considerar trocar por CG com pré-condicionador local

8. **I/O paralelo**
   - Opção simples: cada rank salva `output_rankNNN.vtk`, juntar com ParaView
   - Opção melhor: MPI-IO escrevendo um único `.vtk` paralelo (.pvtk)

#### Arquivos afetados

- [`RTS_main.f90`](sources/RTS_main.f90) — MPI_Init/Finalize
- [`RTS_start.f90`](sources/RTS_start.f90) — topologia, subdomínio
- [`RTS_radiation.f90`](sources/RTS_radiation.f90) — **maior mudança**: sweeps KBA
- [`RTS_solvers.f90`](sources/RTS_solvers.f90) — SOR distribuído + reduções
- [`RTS_bc.f90`](sources/RTS_bc.f90) — diferenciar faces físicas vs internas
- [`RTS_output.f90`](sources/RTS_output.f90) — I/O paralelo

Novo arquivo sugerido:
- `RTS_mpi.f90` — módulo isolando todas as operações MPI (init, halo, sweeps inter-rank)

#### Critério de aceitação

- [ ] Roda com 1 rank (≡ serial) com resultado idêntico
- [ ] Roda com 8 ranks distribuídos (2×2×2) com resultado **dentro da tolerância** do serial
- [ ] Memória por rank cai aproximadamente linearmente com nº de ranks
- [ ] Speedup ≥ **8×** em 16 ranks (em caso 3D grande)
- [ ] Passa em todos os casos de `validation/`

#### Riscos

| Risco | Mitigação |
|-------|-----------|
| KBA é complexo, fácil errar ordens | Implementar octante por octante, validando contra serial |
| Performance ruim por wavefront serial (poucos ranks por face) | Combinar com OpenMP intra-rank na Fase 5 |
| Decadência de convergência do SOR | Trocar por CG/GMRES paralelizado |
| Bug numérico sutil | Comparar campo a campo no benchmark, não só fluxo total |

---

### Fase 5 — Híbrido OpenMP + MPI

**Objetivo:** combinar as duas camadas para extrair máximo dos clusters modernos.
**Esforço:** ~2 semanas
**Pré-requisito:** Fases 2 e 4 funcionando

#### Atividades

1. **Validar coexistência** — compilar com `-fopenmp` + `-lmpi`, rodar com
   `mpirun -n 4 ... OMP_NUM_THREADS=8` (4 ranks × 8 threads = 32 cores)
2. **Tunar layout** — geralmente 1 rank por socket NUMA + threads no socket
3. **Garantir thread-safety nas comunicações MPI** — usar `MPI_THREAD_FUNNELED` (só thread mestre faz MPI) ou superior conforme necessidade
4. **Profiling** — verificar se há sobreposição comunicação/computação (Score-P ou ferramentas similares)

#### Critério de aceitação

- [ ] Speedup total maior que MPI puro ou OpenMP puro isolados
- [ ] Sem deadlocks ou race conditions

#### Riscos

| Risco | Mitigação |
|-------|-----------|
| Conflito entre OpenMP e MPI | Usar `MPI_THREAD_FUNNELED` por simplicidade |
| Carga desbalanceada | Profiling, ajustar mapeamento |

---

### Fase 6 — Validação Completa e Benchmarks

**Objetivo:** documentar o trabalho final com números reais.
**Esforço:** ~2 semanas
**Pré-requisito:** Fase 4 (Fase 5 ideal)

#### Atividades

1. **Rodar todos os casos de `validation/`** e comparar com:
   - Resultado serial original
   - Resultados publicados na tese
2. **Curva de speedup**: tempo vs nº de cores (1, 2, 4, 8, 16, 32, 64, ...)
3. **Curva de eficiência paralela**: speedup / nº de cores
4. **Análise de escalabilidade fraca**: dobra problema E dobra cores → tempo constante?
5. **Análise de memória**: GB por rank em cada configuração
6. **Documentar tudo** em `docs/06-benchmarks.md`

#### Critério de aceitação

- [ ] Todos os benchmarks da tese reproduzidos dentro da tolerância
- [ ] Tabela completa de speedup
- [ ] Identificação clara do "limite de escala" (até onde escala bem)

---

## 5. Ferramentas Necessárias

| Ferramenta | Para que | Fase |
|-----------|----------|------|
| `gfortran` ≥ 10 | Compilação | todas |
| `make` | Build system | todas |
| `mpif90` (OpenMPI ou MPICH) | MPI | 4+ |
| `gprof` ou `perf` | Profiling básico | 0+ |
| `Score-P` ou `Tau` (opcional) | Profiling paralelo avançado | 5+ |
| `ParaView` | Visualizar VTK | 0+ |
| `numpy` + script Python | Comparar campos numéricos | 0+ |
| Acesso a cluster (eventual) | Validar MPI em escala | 4+ |

---

## 6. Caso de Validação Sugerido

**Recomendação:** `validation/3D_Hsu_benchmark/`

**Por quê:**
- 3D (não trivial, exercita os 8 octantes)
- Gas absorvente (exercita WSGG se ativado)
- Tamanho ajustável via `nx,ny,nz` no `input.rts`
- Benchmark conhecido na literatura (Hsu)

**Sugestões de tamanho:**

| Caso | nx=ny=nz | Tempo serial estimado | Memória `IG` serial |
|------|----------|----------------------|---------------------|
| Pequeno (dev) | 32 | segundos | ~70 MB |
| Médio (regressão) | 64 | minutos | ~550 MB |
| Grande (benchmark) | 128 | dezenas de min | ~4 GB |
| Estresse | 256 | horas | ~36 GB |

Começar pequeno (dev rápido), validar com médio, benchmark final no grande.

---

## 7. Métricas de Sucesso

### 7.1 Métricas duras

| Métrica | Meta Fase 1 (OpenMP) | Meta Fase 4 (MPI) |
|---------|---------------------|-------------------|
| Speedup (16 cores) | ≥ 8× | ≥ 12× |
| Speedup (64 cores) | n/a | ≥ 30× |
| Eficiência paralela (16 cores) | ≥ 50% | ≥ 70% |
| Memória por rank (256³, 64 ranks) | n/a | ≤ 2 GB |
| Validação numérica | passa todos | passa todos |

### 7.2 Métricas qualitativas

- [ ] Código continua legível
- [ ] Testes de regressão automatizados
- [ ] Documentação atualizada
- [ ] Comentários nas seções não-óbvias (wavefront, KBA, halo)

---

## 8. Riscos Estratégicos

| Risco | Probabilidade | Impacto | Mitigação |
|-------|--------------|---------|-----------|
| Bug numérico sutil só aparece em escala | Média | Alto | Testes incrementais (1→2→4→8→16 ranks), comparação campo a campo |
| Convergência do SOR degrada | Alta | Médio | Plano B: trocar por CG paralelo |
| Wavefront/KBA escala mal em muitos ranks | Média | Médio | Hibridizar com OpenMP intra-nó |
| Acesso à máquina/cluster atrasa Fase 4 | Média | Alto | Adiantar Fases 1-3 que rodam em desktop |
| Mudanças quebram integração futura com MFSim | Baixa (se manter API) | Alto | Não mexer no nome das rotinas públicas (`rtesolve`, etc.) |
| Pessoa-chave sai do projeto | n/a | Alto | Documentação contínua |

---

## 9. Cronograma Sugerido

```mermaid
gantt
    title Cronograma de Paralelização do RTS
    dateFormat YYYY-MM-DD
    axisFormat %b/%y

    section Fundação
    F0 — Baseline + instrumentação   :f0, 2026-06-01, 1w

    section OpenMP
    F1 — OMP loops fáceis             :f1, after f0, 2w
    F2 — OMP sweeps + red-black       :f2, after f1, 4w

    section Preparação
    F3 — Refactor subdomínio local    :f3, after f2, 3w

    section MPI
    F4 — MPI + KBA                    :f4, after f3, 8w

    section Tuning
    F5 — Híbrido OMP+MPI              :f5, after f4, 2w

    section Validação
    F6 — Benchmarks + docs            :f6, after f5, 2w
```

**Total: ~22 semanas (~5 meses)** assumindo 1 pessoa em dedicação razoável.

> Ajustar conforme o time real. Em paralelo (2 pessoas) pode cair para ~3 meses.

---

## 10. Próximos Passos Imediatos

Para começar **hoje**, na ordem:

1. **Validar acesso ao código:** clonar, abrir, garantir que entendemos a estrutura
2. **Compilar e rodar** o caso `3D_Hsu_benchmark` na configuração atual
3. **Cronometrar** com `cpu_time` quanto demora hoje (será o nosso 0)
4. **Decidir tamanho de benchmark** (sugestão: nx=64 para dev, nx=128 para release)
5. **Criar branch** `paralelizacao` no git para o trabalho
6. **Abrir o documento `05-baseline.md`** (vai capturar a Fase 0)

A partir daí, atacar cada fase em sequência, mantendo as docs em paralelo.

---

## 11. O Que Cada Doc Vai Capturar (futuro)

Conforme as fases avançam, novos docs vão capturar o que foi feito:

| Doc | Conteúdo | Fase associada |
|-----|----------|----------------|
| `05-baseline.md` | Medições atuais (tempo, memória, perfil) | F0 |
| `06-openmp-fase1.md` | Implementação e resultados da OpenMP "fácil" | F1 |
| `07-openmp-fase2.md` | Red-black e wavefront | F2 |
| `08-refactor-local.md` | Refatoração para subdomínio local | F3 |
| `09-mpi-decomposicao.md` | MPI propriamente dito (KBA) | F4 |
| `10-hibrido.md` | OMP + MPI combinados | F5 |
| `11-benchmarks.md` | Curvas finais de speedup, comparação com tese | F6 |

Cada um vai ter: o que mudou no código, antes/depois de performance, problemas encontrados, decisões tomadas.

---

*Documento vivo. Atualizar conforme as fases avançam e decisões são revistas.*
