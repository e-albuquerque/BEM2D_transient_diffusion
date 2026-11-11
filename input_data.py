"""Configuracao dos problemas benchmark para o DCC Moderno (BEM 2D Transiente)."""

# Selecione o problema padrao:
#   1 -> Retangulo 10x8 [-5, 5]x[-4, 4] com u = 1 no contorno (Tese Ex. 2 / Tabela 3.1)
#   2 -> Placa quadrada 6x6 [0, 6]x[0, 6] conducao 1D (Tese Ex. 1 / Matlab dad_0.m)
#   3 -> Disco circular R = 1 com choque termico sucessivo (Tese Ex. 3 / Secao 3.4)
problem = 1


PROBLEMS = {
    1: {
        "name": "retangulo_10x8",
        "file_name": "rect_10x8.geo",
        "title": "Ex. 2 (Tese Sec. 3.3a): Retangulo 10x8 u=1 no contorno [P1]",
        "K_diff": 1.0,
        "tf": 20.0,
        "ninter": 20,  # delta_t = 1.0 s -> coincide com t = 2, 4, ..., 20 s da Tabela 3.1
        "ng": 8,
        "bcs": {
            "bottom": ("temp", 1.0),
            "right": ("temp", 1.0),
            "top": ("temp", 1.0),
            "left": ("temp", 1.0),
        },
        "monitor_points": [(0.0, 0.0)],
        "description": "Retangulo [-5,5]x[-4,4], u_contorno = 1.0, K = 1.0, dt = 1.0 s",
    },
    2: {
        "name": "placa_6x6",
        "file_name": "placa_6x6.geo",
        "title": "Ex. 1 (Tese Sec. 3.2): Placa 6x6 conducao 1D [P2]",
        "K_diff": 1.0,
        "tf": 20.0,
        "ninter": 40,  # delta_t = 0.5 s
        "ng": 8,
        "bcs": {
            "left": ("temp", 300.0),
            "right": ("temp", 0.0),
            "bottom": ("flux", 0.0),
            "top": ("flux", 0.0),
        },
        "monitor_points": [(1.5, 3.0), (3.0, 3.0), (4.5, 3.0)],
        "description": "Placa [0,6]x[0,6], u(0,y)=300, u(6,y)=0, p(x,0)=p(x,6)=0, K=1.0",
    },
    3: {
        "name": "disco_r1",
        "file_name": "disco_r1.geo",
        "title": "Ex. 3 (Tese Sec. 3.4): Disco R=1 choque termico [P3]",
        "K_diff": 1.0,
        "tf": 1.0,
        "ninter": 40,  # delta_t = 0.025 s
        "ng": 8,
        "bcs": {
            # Choque termico duplo: u = 50 para 0 < t <= 0.5 s; u = 100 para t > 0.5 s
            "boundary": ("temp", lambda t: 50.0 if t <= 0.5 + 1.0e-12 else 100.0),
        },
        "monitor_points": [(0.0, 0.0), (0.0, 0.5), (0.0, 0.8)],
        "description": "Disco R=1, u(R,t)=50 (t<=0.5) e 100 (t>0.5), K=1.0, dt=0.025 s",
    },
}


def get_problem_config(problem_id=None):
    """Retorna o dicionario de configuracao para o problema selecionado."""
    pid = problem if problem_id is None else int(problem_id)
    if pid not in PROBLEMS:
        valid = ", ".join(str(k) for k in sorted(PROBLEMS))
        raise ValueError(f"Problema invalido: {pid}. Opcoes validas: {valid}")
    cfg = dict(PROBLEMS[pid])
    cfg["problem_id"] = pid
    cfg["delta_t"] = cfg["tf"] / cfg["ninter"]
    return cfg


# Campos no nivel do modulo para compatibilidade com estilo do projeto de referencia
_default_cfg = get_problem_config(problem)
file_name = _default_cfg["file_name"]
title = _default_cfg["title"]
K_diff = _default_cfg["K_diff"]
tf = _default_cfg["tf"]
ninter = _default_cfg["ninter"]
delta_t = _default_cfg["delta_t"]
ng = _default_cfg["ng"]
bcs = _default_cfg["bcs"]
