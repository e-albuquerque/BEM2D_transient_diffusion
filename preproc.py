"""Pre-processamento com Gmsh para o DCC Moderno (BEM 2D Transiente).

Gera a malha 2D (ordem 1) a partir do arquivo .geo definido em input_data.py,
orienta automaticamente as normais externas dos elementos de contorno usando os
triangulos adjacentes do dominio, extrai pontos internos e escreve dcc_modern.dat.
"""

import argparse
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import meshio
import numpy as np

from input_data import get_problem_config


def generate_gmsh_mesh(geo_path: Path, msh_path: Path) -> None:
    """Executa a API Python do Gmsh para gerar malha 2D linear no formato MSH 2.2."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(geo_path))
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(msh_path))
    finally:
        gmsh.finalize()


def build_edge_to_triangle_map(triangles: np.ndarray) -> dict[tuple[int, int], tuple[int, int, int]]:
    """Mapeia cada aresta (a, b) ordenada para o triangulo (v0, v1, v2) que a contem."""
    edge_map: dict[tuple[int, int], tuple[int, int, int]] = {}
    for tri in triangles:
        v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
        for a, b in ((v0, v1), (v1, v2), (v2, v0)):
            key = (min(a, b), max(a, b))
            edge_map[key] = (v0, v1, v2)
    return edge_map


def orient_boundary_elements(
    points: np.ndarray,
    line_cells: np.ndarray,
    triangles: np.ndarray,
) -> np.ndarray:
    """Garante orientacao consistente para normal externa n = (dy, -dx)/L.

    Em PRECOMPUTE_GEOMETRY (assembly_mod.f90), nx = (y2-y1)/L e ny = -(x2-x1)/L.
    Para que n aponte para fora do dominio, o produto escalar entre n e o vetor
    que vai do ponto medio da aresta ate o vertice oposto do triangulo interno
    adjacente deve ser estritamente negativo.
    """
    edge_map = build_edge_to_triangle_map(triangles)
    oriented = np.empty_like(line_cells)

    for idx, (n1, n2) in enumerate(line_cells):
        a, b = int(n1), int(n2)
        key = (min(a, b), max(a, b))
        if key not in edge_map:
            raise RuntimeError(
                f"Aresta de contorno ({a}, {b}) nao pertence a nenhum triangulo do dominio."
            )
        tri = edge_map[key]
        opp_nodes = [v for v in tri if v not in key]
        if len(opp_nodes) != 1:
            raise RuntimeError(f"Triangulo degenerado {tri} adjacente a aresta {key}.")
        v_opp = opp_nodes[0]

        p1 = points[a]
        p2 = points[b]
        p_opp = points[v_opp]
        mid = 0.5 * (p1 + p2)

        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        # Vetor normal nao unitario segundo PRECOMPUTE_GEOMETRY: (dy, -dx)
        dot = dy * (p_opp[0] - mid[0]) - dx * (p_opp[1] - mid[1])
        if dot > 0.0:
            oriented[idx] = (b, a)
        else:
            oriented[idx] = (a, b)

    return oriented


def extract_boundary_data(mesh: meshio.Mesh, bcs: dict) -> tuple:
    """Extrai elementos de contorno, grupos fisicos e triangulos da malha."""
    points = mesh.points[:, :2].copy()

    line_blocks = []
    line_phys_blocks = []
    tri_blocks = []

    line_phys_all = mesh.cell_data_dict.get("gmsh:physical", {}).get("line", [])
    line_offset = 0

    for cell_block in mesh.cells:
        if cell_block.type == "line":
            n_cells = len(cell_block.data)
            line_blocks.append(cell_block.data.astype(int))
            if len(line_phys_all) >= line_offset + n_cells:
                line_phys_blocks.append(line_phys_all[line_offset : line_offset + n_cells])
            line_offset += n_cells
        elif cell_block.type == "triangle":
            tri_blocks.append(cell_block.data.astype(int))

    if not line_blocks:
        raise RuntimeError("Nenhum elemento linear de contorno ('line') encontrado na malha.")
    if not tri_blocks:
        raise RuntimeError("Nenhum triangulo ('triangle') encontrado na malha.")

    lines = np.vstack(line_blocks)
    triangles = np.vstack(tri_blocks)

    if len(line_phys_blocks) != len(line_blocks):
        raise RuntimeError("Malha sem tags fisicas ('gmsh:physical') para curvas de contorno.")
    line_phys_tags = np.concatenate(line_phys_blocks)

    # Mapeia tag numerica (1D) -> nome fisico do contorno definido no .geo
    tag_to_name = {
        int(tag_info[0]): name
        for name, tag_info in mesh.field_data.items()
        if int(tag_info[1]) == 1
    }

    elem_names = []
    for tag in line_phys_tags:
        name = tag_to_name.get(int(tag))
        if name is None or name not in bcs:
            raise RuntimeError(
                f"Tag fisica {tag} ('{name}') nao encontrada em bcs ({list(bcs.keys())})."
            )
        elem_names.append(name)

    oriented_lines = orient_boundary_elements(points, lines, triangles)
    return points, oriented_lines, triangles, elem_names


def build_compact_numbering(
    points: np.ndarray,
    oriented_lines: np.ndarray,
    triangles: np.ndarray,
) -> tuple:
    """Renumera nos de contorno (1..NN) e extrai nos internos estritos (1..NPI)."""
    boundary_orig_ids = np.unique(oriented_lines.ravel())
    orig_to_boundary_idx = {int(orig): i + 1 for i, orig in enumerate(boundary_orig_ids)}

    boundary_coords = points[boundary_orig_ids]
    conec = np.array(
        [[orig_to_boundary_idx[int(a)], orig_to_boundary_idx[int(b)]] for a, b in oriented_lines],
        dtype=int,
    )

    domain_orig_ids = np.unique(triangles.ravel())
    boundary_set = set(int(i) for i in boundary_orig_ids)
    internal_orig_ids = np.array(
        [int(idx) for idx in domain_orig_ids if int(idx) not in boundary_set],
        dtype=int,
    )
    internal_coords = points[internal_orig_ids] if len(internal_orig_ids) > 0 else np.zeros((0, 2))
    return boundary_coords, conec, internal_coords


def evaluate_bc_values(elem_names: list[str], bcs: dict, delta_t: float, ninter: int) -> tuple:
    """Monta vetor KODE(NE) e matriz VAL(NE, NINTER) ao longo dos passos de tempo."""
    ne = len(elem_names)
    kode = np.zeros(ne, dtype=int)
    vals = np.zeros((ne, ninter), dtype=float)
    times = np.arange(1, ninter + 1, dtype=float) * delta_t

    for i, name in enumerate(elem_names):
        bc_type, bc_spec = bcs[name]
        bc_norm = bc_type.strip().lower()
        if bc_norm in ("temp", "u", "dirichlet", "0"):
            kode[i] = 0
        elif bc_norm in ("flux", "p", "neumann", "1"):
            kode[i] = 1
        else:
            raise ValueError(f"Tipo de contorno desconhecido para '{name}': {bc_type}")

        if callable(bc_spec):
            for k, t in enumerate(times):
                vals[i, k] = float(bc_spec(float(t)))
        else:
            vals[i, :] = float(bc_spec)

    return kode, vals


def write_solver_input(
    dat_path: Path,
    title: str,
    boundary_coords: np.ndarray,
    conec: np.ndarray,
    kode: np.ndarray,
    vals: np.ndarray,
    internal_coords: np.ndarray,
    k_diff: float,
    delta_t: float,
    ng: int,
) -> None:
    """Escreve o arquivo dcc_modern.dat conforme io_mod.f90 (INPUT_DIFF)."""
    nn = len(boundary_coords)
    ne = len(conec)
    ninter = vals.shape[1]
    npi = len(internal_coords)

    with dat_path.open("w", encoding="utf-8") as f:
        f.write(f"{title[:80]}\n")
        f.write(f"{ne:6d} {nn:6d} {ninter:6d} {npi:6d} {k_diff:14.7e} {delta_t:14.7e} {ng:4d}\n")
        for idx, (x, y) in enumerate(boundary_coords, start=1):
            f.write(f"{idx:6d} {x:18.10e} {y:18.10e}\n")
        for idx, ((n1, n2), kd) in enumerate(zip(conec, kode), start=1):
            f.write(f"{idx:6d} {n1:6d} {n2:6d} {kd:4d}\n")
        for row in vals:
            line_vals = " ".join(f"{v:15.7e}" for v in row)
            f.write(f"{line_vals}\n")
        for idx, (x, y) in enumerate(internal_coords, start=1):
            f.write(f"{idx:6d} {x:18.10e} {y:18.10e}\n")


def plot_geometry(
    out_png: Path,
    title: str,
    points: np.ndarray,
    triangles: np.ndarray,
    boundary_coords: np.ndarray,
    conec: np.ndarray,
    kode: np.ndarray,
    internal_coords: np.ndarray,
    monitor_points: list[tuple[float, float]],
) -> None:
    """Plota malha, condicoes de contorno, normais externas e pontos internos."""
    fig, ax = plt.subplots(figsize=(9, 7))

    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    ax.triplot(triang, color="#cccccc", linewidth=0.6, alpha=0.8, zorder=1, label="Malha Gmsh")

    xmin, ymin = boundary_coords.min(axis=0)
    xmax, ymax = boundary_coords.max(axis=0)
    scale = max(xmax - xmin, ymax - ymin)
    normal_len = 0.045 * scale

    dirichlet_labeled = False
    neumann_labeled = False
    normal_x, normal_y, normal_u, normal_v = [], [], [], []

    for (n1, n2), kd in zip(conec, kode):
        p1 = boundary_coords[n1 - 1]
        p2 = boundary_coords[n2 - 1]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        seg_len = np.hypot(dx, dy)
        mid = 0.5 * (p1 + p2)

        normal_x.append(mid[0])
        normal_y.append(mid[1])
        normal_u.append(normal_len * dy / seg_len)
        normal_v.append(-normal_len * dx / seg_len)

        if kd == 0:
            lbl = "Dirichlet ($u$ prescrito)" if not dirichlet_labeled else None
            dirichlet_labeled = True
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#d62728", linewidth=2.4, zorder=3, label=lbl)
        else:
            lbl = "Neumann ($p = \\partial u/\\partial n$ prescrito)" if not neumann_labeled else None
            neumann_labeled = True
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#1f77b4", linewidth=2.4, zorder=3, label=lbl)

    ax.quiver(
        normal_x,
        normal_y,
        normal_u,
        normal_v,
        angles="xy",
        scale_units="xy",
        scale=1.0,
        color="#222222",
        width=0.003,
        headwidth=3.5,
        zorder=4,
        label="Normais externas",
    )

    if len(internal_coords) > 0:
        ax.scatter(
            internal_coords[:, 0],
            internal_coords[:, 1],
            s=14,
            color="#2ca02c",
            alpha=0.85,
            zorder=5,
            label=f"Pontos internos ({len(internal_coords)})",
        )

    if monitor_points:
        mx = [p[0] for p in monitor_points]
        my = [p[1] for p in monitor_points]
        ax.scatter(
            mx,
            my,
            s=90,
            marker="*",
            color="#ff7f0e",
            edgecolors="black",
            linewidths=0.7,
            zorder=6,
            label="Pontos de validacao",
        )

    ax.set_aspect("equal")
    pad = 0.12 * scale
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(f"{title}\nNE = {len(conec)} elementos de contorno | NPI = {len(internal_coords)} pontos internos")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="best", fontsize=8.5, framealpha=0.92)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pre-processamento Gmsh -> dcc_modern.dat")
    parser.add_argument("--problem", type=int, default=None, help="Numero do problema (1, 2 ou 3)")
    args = parser.parse_args()

    work_dir = Path(__file__).resolve().parent
    cfg = get_problem_config(args.problem)

    geo_path = work_dir / cfg["file_name"]
    msh_path = geo_path.with_suffix(".msh")
    dat_path = work_dir / "dcc_modern.dat"
    png_path = work_dir / "geometry.png"

    print(f"[preproc] Problema {cfg['problem_id']}: {cfg['description']}")
    print(f"[preproc] Gerando malha Gmsh a partir de {geo_path.name} -> {msh_path.name}...")
    generate_gmsh_mesh(geo_path, msh_path)

    mesh = meshio.read(msh_path)
    points, oriented_lines, triangles, elem_names = extract_boundary_data(mesh, cfg["bcs"])
    boundary_coords, conec, internal_coords = build_compact_numbering(
        points, oriented_lines, triangles
    )
    kode, vals = evaluate_bc_values(elem_names, cfg["bcs"], cfg["delta_t"], cfg["ninter"])

    write_solver_input(
        dat_path=dat_path,
        title=cfg["title"],
        boundary_coords=boundary_coords,
        conec=conec,
        kode=kode,
        vals=vals,
        internal_coords=internal_coords,
        k_diff=cfg["K_diff"],
        delta_t=cfg["delta_t"],
        ng=cfg["ng"],
    )

    plot_geometry(
        out_png=png_path,
        title=cfg["title"],
        points=points,
        triangles=triangles,
        boundary_coords=boundary_coords,
        conec=conec,
        kode=kode,
        internal_coords=internal_coords,
        monitor_points=cfg["monitor_points"],
    )

    print(
        f"[preproc] Concluido: NE={len(conec)}, NN={len(boundary_coords)}, "
        f"NPI={len(internal_coords)}, NINTER={cfg['ninter']}, dt={cfg['delta_t']:.4g}"
    )
    print(f"[preproc] Arquivos gerados: {msh_path.name}, {dat_path.name}, {png_path.name}")


if __name__ == "__main__":
    main()
