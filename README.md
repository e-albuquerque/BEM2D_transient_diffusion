# DCC Moderno — MEC 2D Transiente para Difusao / Conducao de Calor

Modernizacao do programa Fortran 77 **`DCC`** (tese de mestrado de *Teresinha Costa Effren*, arquivos legados preservados em [`dcc.f`](dcc.f) e [`dcc_orig.f`](dcc_orig.f)) para **Fortran 90 modular**, **alocacao dinamica**, **OpenMP**, solver **LAPACK**, geracao de malha via **Gmsh** e pre/pos-processamento em **Python**.

## Principais Melhorias da Modernizacao

1. **Arquitetura Modular em Fortran 90 (`implicit none` e alocacao dinamica)**:
   - [`data_mod.f90`](data_mod.f90): parametros globais de precisao e vetores/matrizes `allocatable` (sem limites estaticos de 200 nos/passos).
   - [`solver_mod.f90`](solver_mod.f90): interfaces explicitas LAPACK (`DGETRF`, `DGETRS`, `DGESV`) substituindo a eliminacao gaussiana manual (`SLNPD`). A matriz do sistema `C_MAT` e fatorada em $LU$ **uma unica vez** (`FACTOR_SYSTEM`) e reutilizada em todos os `NINTER` passos de tempo (`SOLVE_FACTORED`).
   - [`assembly_mod.f90`](assembly_mod.f90): pre-computacao da geometria e normais do contorno (`PRECOMPUTE_GEOMETRY`), montagem paralela das matrizes transientes (`MATGH`), marcha no tempo (`SOLVE_TRANSIENT`) e avaliacao paralela de pontos internos (`PTINT`) com quadratura de Gauss-Legendre configuravel (`NG`) e transformacao **`sinh`** de Johnston & Elliott (2005) (`SINHTRANS`) para pontos proximos ao contorno.
   - [`io_mod.f90`](io_mod.f90): leitura do formato moderno baseado em conectividade explicita `CONEC(NE, 2)` e exportacao estruturada (`dcc_modern.out`, `ubou.dat`, `pbou.dat`, `uint.dat`).
   - [`main.f90`](main.f90): orquestracao principal com medicao de tempo wall-clock (`omp_get_wtime`) e contagem de threads OpenMP.
2. **Paralelizacao com OpenMP**:
   - Lacos externos de `MATGH` (`collapse(2)` nos pares `(J, I)`) e `PTINT` (pontos internos `KN = 1..NPI`) paralelizados com diretivas OpenMP, mantendo variaveis locais e acumuladores livres de condicoes de corrida.
3. **Precisao Dupla e Correcao da Difusividade $K$**:
   - Compilacao com `-fdefault-real-8`.
   - Correcao da implementacao da formulacao integral transiente (Eq. 2.39–2.40 da tese) na integracao do termo `IQ` para valores gerais de difusividade $K \neq 1$.
4. **Pipeline Python + Gmsh**:
   - [`input_data.py`](input_data.py): selecao dos 3 problemas de validacao (incluindo condicoes de contorno variaveis no tempo).
   - [`preproc.py`](preproc.py): gera malha linear no Gmsh a partir dos arquivos `.geo`, orienta automaticamente as normais externas dos elementos de contorno usando os triangulos do dominio, gera `geometry.png` e grava `dcc_modern.dat`.
   - [`posproc.py`](posproc.py): compara com solucoes analiticas exatas (series de Fourier 1D/2D e funcoes de Bessel $J_0$ com superposicao de Duhamel), calcula erros relativos $L_2$, gera `result_plot.png`, `temperature_map.png` e exporta `resultado_dcc.vtk` para ParaView.

## Requisitos

- **Compilador Fortran + Bibliotecas** (Windows via MSYS2 UCRT64, ou Linux):
  - `gfortran` com suporte OpenMP (`-fopenmp`)
  - `LAPACK` e `BLAS` (`-llapack -lblas`)
- **Python 3** com pacotes:
  - `gmsh`, `meshio`, `numpy`, `matplotlib`

## Compilacao

No Windows (com MSYS2 UCRT64 em `C:\msys64\ucrt64\bin`):
```bat
build.bat
```
Ou via `mingw32-make` / `make`:
```bash
mingw32-make
```

## Fluxo de Execucao

Escolha o problema padrao em [`input_data.py`](input_data.py) (`problem = 1`, `2` ou `3`) ou informe `--problem`:

```bash
python preproc.py --problem 1
dcc_modern.exe
python posproc.py --problem 1
```

## Problemas Benchmark e Resultados de Validacao

| Problema | Geometria (`.geo`) | Malha / Discretizacao | Erro $L_2$ Global | Erro $L_2$ Final | Tempo (8 threads) |
| :--- | :--- | :--- | :---: | :---: | :---: |
| **1** (Tese Ex. 2, Tab. 3.1) | [`rect_10x8.geo`](rect_10x8.geo) `[-5,5]x[-4,4]`, $u=1$ | `NE=68`, `NPI=305`, `dt=1.0 s`, `NINTER=20` | `1.9904%` | `0.0714%` | `0.0480 s` |
| **2** (Tese Ex. 1) | [`placa_6x6.geo`](placa_6x6.geo) `[0,6]x[0,6]`, 1D misto | `NE=56`, `NPI=202`, `dt=0.5 s`, `NINTER=40` | `0.9267%` | `0.0745%` | `0.0550 s` |
| **3** (Tese Ex. 3, Choque) | [`disco_r1.geo`](disco_r1.geo) Disco $R=1$, $u(t)=50 \to 100$ | `NE=48`, `NPI=193`, `dt=0.025 s`, `NINTER=40` | `1.0296%` | `0.0942%` | `0.0480 s` |

Em todos os casos, `DCC Moderno` com 8 threads OpenMP reduz o tempo no Problema 2 de `0.2420 s` (1 thread) para `0.0570 s` (aceleracao ~4.2x), mantendo precisao consistente tambem para $K \neq 1$ (validado em `K = 2.5`, `dt = 0.2 s` com mesmo erro global `0.9267%`).

## Saidas Visuais Geradas

- `geometry.png`: malha triangular do Gmsh, contorno Dirichlet/Neumann, normais externas e pontos internos.
- `result_plot.png`: historico temporal nos pontos de validacao + perfis espaciais (inclui comparacao com a Tabela 3.1 no Problema 1).
- `temperature_map.png`: mapas de contorno 2D da evolucao da temperatura em 4 instantes de tempo.
- `resultado_dcc.vtk`: malha e campos de temperatura transiente para visualizacao no ParaView.
