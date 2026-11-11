"""Pos-processamento, validacao analitica e visualizacao do DCC Moderno.

Le os resultados gerados por dcc_modern.exe (ubou.dat, pbou.dat, uint.dat),
compara com solucoes analiticas exatas (e Tabela 3.1 da tese quando aplicavel),
imprime erros relativos L2, gera graficos de validacao (result_plot.png e
temperature_map.png) e exporta campo transiente para ParaView (resultado_dcc.vtk).
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import meshio
import numpy as np

from input_data import get_problem_config


# ---------------------------------------------------------------------------
# Funcoes especiais (Bessel J0, J1 e zeros de J0 em NumPy puro)
# ---------------------------------------------------------------------------

def _bessel_j0_j1(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calcula J0(x) e J1(x) com precisao de maquina via quadratura periodica."""
    x_arr = np.asarray(x, dtype=float)
    theta = np.linspace(0.0, 2.0 * np.pi, 1024, endpoint=False)
    arg = x_arr[..., None] * np.cos(theta)
    j0 = np.mean(np.cos(arg), axis=-1)
    j1 = np.mean(np.sin(arg) * np.cos(theta), axis=-1)
    return j0, j1


def _j0_zeros(n_roots: int = 60) -> np.ndarray:
    """Encontra os primeiros n_roots zeros positivos de J0(alpha) = 0 via Newton."""
    idx = np.arange(1, n_roots + 1, dtype=float)
    beta = (idx - 0.25) * np.pi
    # Aproximacao assintotica de McMahon
    alpha = beta + 1.0 / (8.0 * beta) - 31.0 / (384.0 * beta**3)
    for _ in range(7):
        j0_val, j1_val = _bessel_j0_j1(alpha)
        # J0'(alpha) = -J1(alpha) => passo de Newton: +J0/J1
        alpha = alpha + j0_val / j1_val
    return alpha


# ---------------------------------------------------------------------------
# Solucoes analiticas exatas dos problemas benchmark
# ---------------------------------------------------------------------------

def analytical_rect_10x8(
    x: np.ndarray, y: np.ndarray, t: float, k_diff: float = 1.0, u0: float = 1.0
) -> np.ndarray:
    """Solucao analitica 2D por serie de Fourier no retangulo [-5,5]x[-4,4]."""
    if t <= 0.0:
        return np.zeros_like(x, dtype=float)
    a, b = 5.0, 4.0
    modes = np.arange(0, 60, dtype=float)
    odd = 2.0 * modes + 1.0
    sign = (-1.0) ** modes

    kx = odd * np.pi / (2.0 * a)
    ky = odd * np.pi / (2.0 * b)

    # Termos separaveis em x e y com pesos espaciais
    wx = (sign / odd)[:, None] * np.cos(kx[:, None] * x[None, :])
    wy = (sign / odd)[:, None] * np.cos(ky[:, None] * y[None, :])
    ex = np.exp(-k_diff * (kx**2) * t)[:, None]
    ey = np.exp(-k_diff * (ky**2) * t)[:, None]

    sx = np.sum(wx * ex, axis=0)
    sy = np.sum(wy * ey, axis=0)
    return u0 * (1.0 - (16.0 / (np.pi**2)) * sx * sy)


def analytical_placa_6x6(
    x: np.ndarray, _y: np.ndarray, t: float, k_diff: float = 1.0, u0: float = 300.0
) -> np.ndarray:
    """Solucao analitica 1D por serie de Fourier na placa [0,6]x[0,6] (Tese Eq. A2.19)."""
    if t <= 0.0:
        return np.zeros_like(x, dtype=float)
    length = 6.0
    n = np.arange(1, 151, dtype=float)[:, None]
    x_row = np.asarray(x, dtype=float)[None, :]
    series = np.sum(
        (1.0 / n)
        * np.sin(n * np.pi * x_row / length)
        * np.exp(-k_diff * (n * np.pi / length) ** 2 * t),
        axis=0,
    )
    return u0 * (1.0 - np.asarray(x, dtype=float) / length) - (2.0 * u0 / np.pi) * series


_DISK_ALPHAS = None
_DISK_J1 = None


def _analytical_disk_step(
    r: np.ndarray, tau: float, k_diff: float, u_step: float, radius: float = 1.0
) -> np.ndarray:
    """Resposta a degrau unitario de temperatura no contorno de um disco de raio R."""
    global _DISK_ALPHAS, _DISK_J1
    if tau <= 0.0:
        return np.zeros_like(r, dtype=float)
    if _DISK_ALPHAS is None:
        _DISK_ALPHAS = _j0_zeros(60)
        _, _DISK_J1 = _bessel_j0_j1(_DISK_ALPHAS)

    alphas = _DISK_ALPHAS
    j1_vals = _DISK_J1
    r_arr = np.clip(np.asarray(r, dtype=float), 0.0, radius)

    arg = alphas[:, None] * (r_arr[None, :] / radius)
    j0_vals, _ = _bessel_j0_j1(arg)
    decay = np.exp(-k_diff * (alphas**2) * tau / (radius**2))[:, None]
    coeff = (1.0 / (alphas * j1_vals))[:, None]
    return u_step * (1.0 - 2.0 * np.sum(coeff * j0_vals * decay, axis=0))


def analytical_disco_r1(
    x: np.ndarray, y: np.ndarray, t: float, k_diff: float = 1.0
) -> np.ndarray:
    """Solucao analitica por superposicao de Duhamel no disco R=1 (Tese Sec. 3.4)."""
    r = np.hypot(x, y)
    u = _analytical_disk_step(r, t, k_diff, 50.0, radius=1.0)
    if t > 0.5 + 1.0e-12:
        u = u + _analytical_disk_step(r, t - 0.5, k_diff, 50.0, radius=1.0)
    return u


ANALYTICAL_SOLVERS = {
    1: analytical_rect_10x8,
    2: analytical_placa_6x6,
    3: analytical_disco_r1,
}

# Dados da Tabela 3.1 (Tese Teresinha Effren, pag. 50) para o centro A(0,0) do Ex. 2
THESIS_TABLE_3_1 = {
    "t": np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float),
    "anal": np.array([0.114, 0.420, 0.646, 0.786, 0.871, 0.922, 0.953, 0.972, 0.983, 0.990]),
    "dcc_28": np.array([0.092, 0.388, 0.621, 0.769, 0.859, 0.915, 0.948, 0.968, 0.981, 0.988]),
}


# ---------------------------------------------------------------------------
# Leitura de entrada, saidas do solver e reconstrucao de campos nos nos
# ---------------------------------------------------------------------------

def read_solver_input(dat_path: Path) -> dict:
    """Le cabecalho e discretizacao do arquivo dcc_modern.dat."""
    with dat_path.open("r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    title = lines[0]
    parts = lines[1].split()
    ne = int(parts[0])
    nn = int(parts[1])
    ninter = int(parts[2])
    npi = int(parts[3])
    k_diff = float(parts[4])
    delta_t = float(parts[5])
    ng = int(parts[6])

    offset = 2
    b_coords = np.array(
        [[float(p) for p in lines[offset + i].split()[1:3]] for i in range(nn)],
        dtype=float,
    )
    offset += nn

    conec = np.array(
        [[int(p) for p in lines[offset + i].split()[1:3]] for i in range(ne)],
        dtype=int,
    )
    kode = np.array(
        [int(lines[offset + i].split()[3]) for i in range(ne)],
        dtype=int,
    )
    return {
        "title": title,
        "ne": ne,
        "nn": nn,
        "ninter": ninter,
        "npi": npi,
        "k_diff": k_diff,
        "delta_t": delta_t,
        "ng": ng,
        "boundary_coords": b_coords,
        "conec": conec,
        "kode": kode,
    }


def detect_problem_id(cli_problem: int | None, job_title: str) -> int:
    """Determina o problema via argumento CLI ou marcador [P#] no titulo."""
    if cli_problem is not None:
        return int(cli_problem)
    match = re.search(r"\[P(\d+)\]", job_title)
    if match:
        return int(match.group(1))
    return int(get_problem_config()["problem_id"])


def read_fortran_metadata(out_path: Path) -> dict[str, str]:
    """Extrai contagem de threads OpenMP e tempo de execucao de dcc_modern.out."""
    meta = {}
    if not out_path.exists():
        return meta
    text = out_path.read_text(encoding="utf-8", errors="ignore")
    m_thr = re.search(r"Threads OpenMP:\s*(\d+)", text)
    m_sec = re.search(r"Tempo total \(wall-clock\):\s*([\d.E+-]+)\s*s", text)
    if m_thr:
        meta["threads"] = m_thr.group(1)
    if m_sec:
        meta["elapsed"] = m_sec.group(1)
    return meta


def reconstruct_nodal_temperatures(
    mesh: meshio.Mesh,
    ubou_vals: np.ndarray,
    kode: np.ndarray,
    uint_coords: np.ndarray,
    uint_vals: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Projeta u_elem (contorno) e u_int (interno) para todos os vertices da malha."""
    points = mesh.points[:, :2].copy()
    line_blocks = [cb.data.astype(int) for cb in mesh.cells if cb.type == "line"]
    tri_blocks = [cb.data.astype(int) for cb in mesh.cells if cb.type == "triangle"]
    lines = np.vstack(line_blocks)
    triangles = np.vstack(tri_blocks)

    boundary_orig_ids = np.unique(lines.ravel())
    boundary_set = set(int(i) for i in boundary_orig_ids)
    domain_orig_ids = np.unique(triangles.ravel())
    internal_orig_ids = np.array(
        [int(idx) for idx in domain_orig_ids if int(idx) not in boundary_set],
        dtype=int,
    )

    if len(internal_orig_ids) != len(uint_coords):
        raise RuntimeError(
            f"Numero de nos internos diverge: malha={len(internal_orig_ids)} vs uint.dat={len(uint_coords)}"
        )

    # Nos vertices de contorno com mistura Dirichlet/Neumann, prioriza elementos Dirichlet
    n_pts = len(points)
    ninter = ubou_vals.shape[1]
    u_nodes = np.zeros((n_pts, ninter), dtype=float)
    d_accum = np.zeros((n_pts, ninter), dtype=float)
    d_counts = np.zeros(n_pts, dtype=float)
    b_accum = np.zeros((n_pts, ninter), dtype=float)
    b_counts = np.zeros(n_pts, dtype=float)

    for e_idx, (n1, n2) in enumerate(lines):
        if kode[e_idx] == 0:
            d_accum[n1] += ubou_vals[e_idx]
            d_accum[n2] += ubou_vals[e_idx]
            d_counts[n1] += 1.0
            d_counts[n2] += 1.0
        b_accum[n1] += ubou_vals[e_idx]
        b_accum[n2] += ubou_vals[e_idx]
        b_counts[n1] += 1.0
        b_counts[n2] += 1.0

    mask_d = d_counts > 0
    mask_n = (b_counts > 0) & (~mask_d)
    u_nodes[mask_d] = d_accum[mask_d] / d_counts[mask_d, None]
    u_nodes[mask_n] = b_accum[mask_n] / b_counts[mask_n, None]
    u_nodes[internal_orig_ids] = uint_vals

    return points, triangles, u_nodes


def find_nearest_internal_point(xi: np.ndarray, yi: np.ndarray, target: tuple[float, float]) -> int:
    """Encontra o indice do ponto interno mais proximo de target=(xt, yt)."""
    dists = np.hypot(xi - target[0], yi - target[1])
    return int(np.argmin(dists))


# ---------------------------------------------------------------------------
# Graficos e VTK
# ---------------------------------------------------------------------------

def plot_validation(
    out_png: Path,
    problem_id: int,
    cfg: dict,
    times: np.ndarray,
    xi: np.ndarray,
    yi: np.ndarray,
    uint_vals: np.ndarray,
    points: np.ndarray,
    triangles: np.ndarray,
    u_nodes: np.ndarray,
    k_diff: float,
) -> None:
    """Gera result_plot.png com comparacao temporal e perfil espacial BEM vs Analitico."""
    anal_fn = ANALYTICAL_SOLVERS[problem_id]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.4))

    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    markers = ["o", "s", "^"]

    # --- Painel esquerdo: serie temporal nos pontos de monitoramento ---
    t_dense = np.linspace(0.0, times[-1], 250)
    for idx, pt in enumerate(cfg["monitor_points"]):
        p_idx = find_nearest_internal_point(xi, yi, pt)
        px, py = xi[p_idx], yi[p_idx]
        c = colors[idx % len(colors)]
        m = markers[idx % len(markers)]

        u_exact_dense = np.array([anal_fn(np.array([px]), np.array([py]), float(t), k_diff)[0] for t in t_dense])
        label_base = f"({px:.2g}, {py:.2g})"

        ax1.plot(t_dense, u_exact_dense, color=c, linewidth=1.8, label=f"Analitico {label_base}")
        step_skip = 1 if len(times) <= 25 else 2
        ax1.plot(
            times[::step_skip],
            uint_vals[p_idx, ::step_skip],
            linestyle="none",
            marker=m,
            markersize=5.2,
            color=c,
            markerfacecolor="white",
            markeredgewidth=1.4,
            label=f"DCC Moderno {label_base}",
        )

    if problem_id == 1:
        ax1.plot(
            THESIS_TABLE_3_1["t"],
            THESIS_TABLE_3_1["dcc_28"],
            linestyle="none",
            marker="x",
            markersize=6.5,
            color="#444444",
            markeredgewidth=1.5,
            label="DCC Tese (28 elem., Tabela 3.1)",
        )

    ax1.set_xlabel("Tempo $t$ (s)")
    ax1.set_ylabel("Temperatura $u(x, y, t)$")
    ax1.set_title("Historico Temporal nos Pontos de Validacao")
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="best", fontsize=8.0, framealpha=0.92)

    # --- Painel direito: perfil espacial em instantes selecionados ---
    ninter = len(times)
    if problem_id == 1:
        snap_indices = [1, 5, 11, ninter - 1]
        triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
        x_line = np.linspace(-5.0, 5.0, 120)
        y_line = np.zeros_like(x_line)
        x_samp = np.linspace(-4.7, 4.7, 17)
        y_samp = np.zeros_like(x_samp)
        for s_i, k_step in enumerate(snap_indices):
            c = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"][s_i]
            t_val = float(times[k_step])
            u_ex = anal_fn(x_line, y_line, t_val, k_diff)
            interp = mtri.LinearTriInterpolator(triang, u_nodes[:, k_step])
            u_bem = np.asarray(interp(x_samp, y_samp), dtype=float)
            ax2.plot(x_line, u_ex, color=c, linewidth=1.8, label=f"$t={t_val:.0f}$ s (Analitico)")
            ax2.plot(x_samp, u_bem, "o", color=c, markersize=4.5, markerfacecolor="white", markeredgewidth=1.3)
        ax2.set_xlabel("Coordenada $x$ em $y = 0$")
        ax2.set_ylabel("Temperatura $u(x, 0, t)$")
        ax2.set_title("Perfil Espacial ao Longo de $y = 0$ (simbolos: BEM)")

    elif problem_id == 2:
        snap_indices = [3, 9, 19, ninter - 1]
        x_line = np.linspace(0.0, 6.0, 150)
        y_line = np.zeros_like(x_line)
        triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
        mask_mid = np.abs(yi - 3.0) < 0.7
        for s_i, k_step in enumerate(snap_indices):
            c = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"][s_i]
            t_val = float(times[k_step])
            u_ex = anal_fn(x_line, y_line, t_val, k_diff)
            ax2.plot(x_line, u_ex, color=c, linewidth=1.8, label=f"$t={t_val:.1f}$ s")
            ax2.plot(
                xi[mask_mid],
                uint_vals[mask_mid, k_step],
                "o",
                color=c,
                markersize=4.2,
                markerfacecolor="white",
                markeredgewidth=1.2,
            )
        ax2.set_xlabel("Coordenada $x$")
        ax2.set_ylabel("Temperatura $u(x, t)$")
        ax2.set_title("Perfil 1D ao Longo de $x$ (simbolos: BEM)")

    else:
        snap_indices = [3, 11, 19, 27, ninter - 1]
        r_vals = np.hypot(xi, yi)
        r_grid = np.linspace(0.0, 1.0, 150)
        palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#d62728"]
        for s_i, k_step in enumerate(snap_indices):
            c = palette[s_i]
            t_val = float(times[k_step])
            u_ex = anal_fn(r_grid, np.zeros_like(r_grid), t_val, k_diff)
            ax2.plot(r_grid, u_ex, color=c, linewidth=1.8, label=f"$t={t_val:.3f}$ s")
            ax2.plot(
                r_vals[::2],
                uint_vals[::2, k_step],
                "o",
                color=c,
                markersize=3.8,
                alpha=0.85,
                markerfacecolor="white",
                markeredgewidth=1.0,
            )
        ax2.set_xlabel("Raio $r = \\sqrt{x^2 + y^2}$")
        ax2.set_ylabel("Temperatura $u(r, t)$")
        ax2.set_title("Perfil Radial nos Pontos Internos (simbolos: BEM)")

    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="best", fontsize=8.0, framealpha=0.92)

    fig.suptitle(f"Validacao DCC Moderno — {cfg['description']}", fontsize=11.5, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_temperature_maps(
    out_png: Path,
    cfg: dict,
    times: np.ndarray,
    points: np.ndarray,
    triangles: np.ndarray,
    u_nodes: np.ndarray,
) -> None:
    """Gera mapa de contorno 2D de temperatura em 4 instantes de tempo."""
    ninter = len(times)
    snap_indices = sorted(
        {
            max(0, int(round(0.08 * (ninter - 1)))),
            max(1, int(round(0.30 * (ninter - 1)))),
            max(2, int(round(0.60 * (ninter - 1)))),
            ninter - 1,
        }
    )[:4]

    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    vmin = float(np.min(u_nodes[:, snap_indices]))
    vmax = float(np.max(u_nodes[:, snap_indices]))
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0
    levels = np.linspace(vmin, vmax, 25)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0), constrained_layout=True)
    cf = None
    for ax, k_step in zip(axes.ravel(), snap_indices):
        field = u_nodes[:, k_step]
        cf = ax.tricontourf(triang, field, levels=levels, cmap="inferno", extend="neither")
        ax.tricontour(triang, field, levels=10, colors="white", linewidths=0.45, alpha=0.55)
        ax.triplot(triang, color="#000000", linewidth=0.22, alpha=0.18)
        ax.set_aspect("equal")
        ax.set_title(f"$t = {times[k_step]:.3g}$ s  (passo {k_step + 1}/{ninter})", fontsize=10)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")

    cbar = fig.colorbar(cf, ax=axes, orientation="vertical", shrink=0.88, pad=0.02)
    cbar.set_label("Temperatura $u(x, y, t)$")
    fig.suptitle(f"Evolucao do Campo de Temperatura 2D — {cfg['title']}", fontsize=12, fontweight="bold")
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def export_vtk(
    out_vtk: Path,
    times: np.ndarray,
    points: np.ndarray,
    triangles: np.ndarray,
    u_nodes: np.ndarray,
    anal_final: np.ndarray,
) -> None:
    """Exporta snapshots temporais de temperatura e erro final em formato VTK."""
    ninter = len(times)
    n_snaps = min(10, ninter)
    snap_indices = np.unique(np.linspace(0, ninter - 1, n_snaps, dtype=int))

    pts_3d = np.column_stack([points, np.zeros(len(points))])
    point_data = {}
    for k_step in snap_indices:
        label = f"T_step_{k_step + 1:03d}_t_{times[k_step]:.3f}s".replace(".", "p")
        point_data[label] = u_nodes[:, k_step]
    point_data["T_final"] = u_nodes[:, -1]
    point_data["abs_error_final"] = np.abs(u_nodes[:, -1] - anal_final)

    vtk_mesh = meshio.Mesh(
        points=pts_3d,
        cells=[("triangle", triangles)],
        point_data=point_data,
    )
    meshio.write(str(out_vtk), vtk_mesh, file_format="vtk")


# ---------------------------------------------------------------------------
# Rotina principal
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Pos-processamento e validacao do DCC Moderno")
    parser.add_argument("--problem", type=int, default=None, help="Numero do problema (1, 2 ou 3)")
    args = parser.parse_args()

    work_dir = Path(__file__).resolve().parent
    dat_info = read_solver_input(work_dir / "dcc_modern.dat")
    problem_id = detect_problem_id(args.problem, dat_info["title"])
    cfg = get_problem_config(problem_id)

    msh_path = (work_dir / cfg["file_name"]).with_suffix(".msh")
    mesh = meshio.read(msh_path)

    ubou_raw = np.loadtxt(work_dir / "ubou.dat", dtype=float)
    uint_raw = np.loadtxt(work_dir / "uint.dat", dtype=float)
    if ubou_raw.ndim == 1:
        ubou_raw = ubou_raw[None, :]
    if uint_raw.ndim == 1:
        uint_raw = uint_raw[None, :]

    ubou_vals = ubou_raw[:, 2:]
    xi, yi = uint_raw[:, 0], uint_raw[:, 1]
    uint_vals = uint_raw[:, 2:]

    ninter = dat_info["ninter"]
    delta_t = dat_info["delta_t"]
    k_diff = dat_info["k_diff"]
    times = np.arange(1, ninter + 1, dtype=float) * delta_t

    points, triangles, u_nodes = reconstruct_nodal_temperatures(
        mesh, ubou_vals, dat_info["kode"], uint_raw[:, :2], uint_vals
    )

    anal_fn = ANALYTICAL_SOLVERS[problem_id]
    u_anal_int = np.column_stack([anal_fn(xi, yi, float(t), k_diff) for t in times])
    anal_final_nodes = anal_fn(points[:, 0], points[:, 1], float(times[-1]), k_diff)

    global_l2 = np.linalg.norm(uint_vals - u_anal_int) / np.linalg.norm(u_anal_int)
    final_l2 = np.linalg.norm(uint_vals[:, -1] - u_anal_int[:, -1]) / np.linalg.norm(u_anal_int[:, -1])
    max_abs = float(np.max(np.abs(uint_vals - u_anal_int)))

    meta = read_fortran_metadata(work_dir / "dcc_modern.out")
    print(f"[posproc] Problema {problem_id}: {cfg['title']}")
    if meta:
        print(
            f"[posproc] Execucao Fortran: {meta.get('threads', '?')} threads OpenMP, "
            f"tempo = {meta.get('elapsed', '?')} s"
        )
    print(f"[posproc] Erro L2 relativo global (todos os nos internos x passos): {global_l2 * 100:.4f}%")
    print(f"[posproc] Erro L2 relativo no instante final (t={times[-1]:.3g} s): {final_l2 * 100:.4f}%")
    print(f"[posproc] Erro absoluto maximo interno: {max_abs:.6f}")

    if problem_id == 1:
        p0 = find_nearest_internal_point(xi, yi, (0.0, 0.0))
        step_idx = np.array([1, 3, 5, 7, 9, 11, 13, 15, 17, 19], dtype=int)
        u_bem_tab = uint_vals[p0, step_idx]
        max_tab_diff = float(np.max(np.abs(u_bem_tab - THESIS_TABLE_3_1["anal"])))
        print(f"[posproc] Centro A(0,0) vs Tabela 3.1 (analitico) max diff = {max_tab_diff:.6f}")
        print("   t(s)   Anal.(Tese)   DCC(Tese 28)   DCC Moderno")
        for ti, ua, ut, um in zip(
            THESIS_TABLE_3_1["t"],
            THESIS_TABLE_3_1["anal"],
            THESIS_TABLE_3_1["dcc_28"],
            u_bem_tab,
        ):
            print(f"  {ti:4.0f}     {ua:8.3f}      {ut:8.3f}      {um:8.4f}")

    plot_validation(
        out_png=work_dir / "result_plot.png",
        problem_id=problem_id,
        cfg=cfg,
        times=times,
        xi=xi,
        yi=yi,
        uint_vals=uint_vals,
        points=points,
        triangles=triangles,
        u_nodes=u_nodes,
        k_diff=k_diff,
    )
    plot_temperature_maps(
        out_png=work_dir / "temperature_map.png",
        cfg=cfg,
        times=times,
        points=points,
        triangles=triangles,
        u_nodes=u_nodes,
    )
    export_vtk(
        out_vtk=work_dir / "resultado_dcc.vtk",
        times=times,
        points=points,
        triangles=triangles,
        u_nodes=u_nodes,
        anal_final=anal_final_nodes,
    )
    print("[posproc] Gerados: result_plot.png, temperature_map.png, resultado_dcc.vtk")


if __name__ == "__main__":
    main()
