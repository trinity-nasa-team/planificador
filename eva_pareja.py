"""eva_pareja.py - EVA en pareja y radio con el hábitat (sobre Marswalk v4, sin cambiar mision_v3.py).

Dos cosas que revisa a lo largo de una EVA ya trazada (hábitat -> paradas -> hábitat):

1. PAREJA (en la NASA nadie sale solo)
   - Oxígeno de cada astronauta, con su propio peso y su propio tanque (modelo metabólico de mision_v3:
     Pandolf + Santee, traje de 104 kg, gravedad de Marte). Punto de no retorno = primer punto desde el que a
     alguno de los dos ya no le alcanza el O2 para volver (guardando la reserva).
   - RESCATE ("compañero lesionado"): en cada punto se supone que uno de los dos se lesiona y el otro lo ayuda a
     volver. Se prueban los dos papeles y se toma el peor. Ambos deben llegar al hábitat con O2 (aquí sí se puede
     usar la reserva) y dentro de las horas del traje. `sin_rescate` = primer punto desde el que ya no alcanzan:
     la misión debería cambiarse por una más cercana o más corta.
     SUPUESTOS PROPIOS (no son datos de la NASA): con un lesionado se camina a VEL_ASISTIDA de la velocidad normal
     y quien ayuda carga CARGA_AYUDA_KG extra (parte del peso del compañero).

2. RADIO
   - Línea de vista entre la antena del hábitat (ANTENA_HABITAT_M) y la del traje (ANTENA_TRAJE_M) sobre el
     relieve, con la curvatura de Marte. Con 10 m y 2 m el horizonte queda a ~12 km: dentro de unos km lo que corta
     la señal son las lomas y los escarpes.
   - con_radio(): copia del mapa en la que las celdas sin señal cuestan W_SIN_SENAL más, para que las rutas las
     eviten cuando pueden.

Uso:
    import eva_pareja as ep
    pareja = [{"masa_kg": 80, "o2_L": 380}, {"masa_kg": 70, "o2_L": 380}]
    r = ep.perfil_pareja(mapa_astronauta, camino, paradas, habitat, pareja)
    cob = ep.cobertura_radio(mapa_astronauta, habitat, radio_m=4000)
    ep.senal_en_camino(camino, cob)
"""
import math
from dataclasses import replace

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

import mision_v3 as v3

R_MARTE = 3396190.0
O2_TRAJE_L = 380.0          # tanque primario del traje EMU de la NASA: 0.54 kg de O2 ~ 380 L
HORAS_TRAJE = 8.0           # soporte vital primario del EMU
RESERVA = 0.10              # en la EVA normal se guarda 10 % del tanque; en un rescate se puede usar

VEL_ASISTIDA = 0.5          # SUPUESTO: con un compañero lesionado se camina a la mitad de la velocidad
CARGA_AYUDA_KG = 30.0       # SUPUESTO: quien ayuda soporta ~30 kg-masa del compañero

ANTENA_HABITAT_M = 10.0     # mástil del hábitat
ANTENA_TRAJE_M = 2.0        # antena en la mochila del traje
W_SIN_SENAL = 3.0           # costo extra por celda sin señal (como una pendiente de ~25° en el mapa de costos)

PAREJA_DEFECTO = ({"masa_kg": 80.0, "o2_L": O2_TRAJE_L}, {"masa_kg": 70.0, "o2_L": O2_TRAJE_L})


# ---------------------------------------------------------------- pareja
def astronauta(masa_kg=80.0, kit_kg=None):
    """mision_v3.ASTRONAUTA con otra masa corporal (y, si se da, otro kit)."""
    cambios = {"masa_kg": float(masa_kg)}
    if kit_kg is not None:
        cambios["kit_kg"] = float(kit_kg)
    return replace(v3.ASTRONAUTA, **cambios)


def _aristas(mapa, base, radio_m):
    """Grafo de 8 vecinos alrededor de la base (mismas reglas que costo_regreso), con lo necesario para pesarlo."""
    v = v3.Ventana.alrededor(mapa, [base], margen_m=radio_m, centro=base, radio_m=radio_m)
    n, m = v.costo.shape
    libre = np.isfinite(v.costo)
    dem = np.nan_to_num(v.sl(mapa.dem))
    arena = v.sl(mapa.arena) if mapa.arena is not None else np.zeros((n, m), bool)
    lento = v.sl(mapa.riesgo | mapa.duna)
    filas, cols = np.meshgrid(np.arange(n), np.arange(m), indexing="ij")
    A = {k: [] for k in ("o", "d", "dh", "dz", "lento", "arena")}
    for df, dc in v3.VECINOS:
        sl_a = (slice(max(0, -df), min(n, n - df)), slice(max(0, -dc), min(m, m - dc)))
        fa, ca = filas[sl_a], cols[sl_a]
        fb, cb = fa + df, ca + dc
        ok = libre[fa, ca] & libre[fb, cb]
        if df and dc:
            ok &= libre[fa, cb] & libre[fb, ca]
        fa, ca, fb, cb = fa[ok], ca[ok], fb[ok], cb[ok]
        A["o"].append(fa * m + ca)
        A["d"].append(fb * m + cb)
        A["dh"].append(np.full(len(fa), v3.RES * math.hypot(df, dc)))
        A["dz"].append(dem[fb, cb] - dem[fa, ca])
        A["lento"].append(lento[fb, cb])
        A["arena"].append(arena[fb, cb])
    return v, (n, m), {k: np.concatenate(x) for k, x in A.items()}


def regresos(mapa, base, agentes, radio_m, vel_factor=1.0):
    """Energía metabólica mínima [Wh] para volver a `base` desde cada celda, para varios agentes a la vez
    (`agentes` = {clave: Agente}), y el tiempo mínimo [h] en la clave "h". vel_factor < 1 = regreso más lento.
    Arreglos del tamaño del mapa (NaN = sin regreso o fuera del radio). Con un solo agente y vel_factor = 1
    coincide con mision_v3.costo_regreso."""
    v, (n, m), A = _aristas(mapa, base, radio_m)
    i0 = v.local(base)

    def dijk(w):
        g = csr_matrix((w, (A["d"], A["o"])), shape=(n * m, n * m))      # aristas invertidas: desde la base
        dist = dijkstra(g, directed=True, indices=i0)
        arr = np.full(mapa.forma, np.nan)
        v.sl(arr)[...] = np.where(np.isfinite(dist), dist, np.nan).reshape(n, m)
        return arr

    salida = {}
    for clave, ag in agentes.items():
        vel = np.where(A["lento"], ag.vel_riesgo_m_h, ag.vel_seguro_m_h) * vel_factor
        t_h = np.hypot(A["dh"], A["dz"]) / vel
        eta = np.where(A["arena"], ag.eta_arena, ag.eta_terreno)
        salida[clave] = dijk(v3.pandolf_w(ag, vel / 3600, 100 * A["dz"] / A["dh"], eta) * t_h)
        if "h" not in salida:
            salida["h"] = dijk(t_h)
    return salida


def tiempos_camino(mapa, camino, agente, paradas=()):
    """Horas transcurridas desde la salida en cada punto del camino (con las paradas de estudio)."""
    a = np.asarray(camino)
    dz = np.diff(np.nan_to_num(mapa.dem[a[:, 0], a[:, 1]]))
    dh = v3.RES * np.hypot(*np.diff(a, axis=0).T)
    lento = (mapa.riesgo | mapa.duna)[a[1:, 0], a[1:, 1]]
    vel = np.where(lento, agente.vel_riesgo_m_h, agente.vel_seguro_m_h)
    t = np.r_[0.0, np.cumsum(np.hypot(dh, dz) / vel)]
    for i in paradas:
        t[i:] += agente.horas_por_objetivo
    return t


def distancias_camino(mapa, camino):
    """Metros recorridos (3D) hasta cada punto del camino."""
    a = np.asarray(camino)
    dz = np.diff(np.nan_to_num(mapa.dem[a[:, 0], a[:, 1]]))
    return np.r_[0.0, np.cumsum(np.hypot(v3.RES * np.hypot(*np.diff(a, axis=0).T), dz))]


def regresos_pareja(mapa, base, pareja, radio_m):
    """Los dos juegos de regresos que necesita perfil_pareja: normal (cada uno a su paso) y asistido."""
    ag = [astronauta(p.get("masa_kg", 80.0), p.get("kit_kg")) for p in pareja]
    normal = regresos(mapa, base, {0: ag[0], 1: ag[1]}, radio_m)
    asistido = regresos(mapa, base, {
        "ayuda0": replace(ag[0], kit_kg=ag[0].kit_kg + CARGA_AYUDA_KG),
        "ayuda1": replace(ag[1], kit_kg=ag[1].kit_kg + CARGA_AYUDA_KG),
        "herido0": ag[0], "herido1": ag[1]}, radio_m, VEL_ASISTIDA)
    return normal, asistido


def perfil_pareja(mapa, camino, paradas, base, pareja=PAREJA_DEFECTO, horas_max=HORAS_TRAJE, reserva=RESERVA,
                  radio_m=None, regresos_calculados=None):
    """Revisa una EVA ya trazada (`camino` en celdas, empieza y termina en `base`; `paradas` = índices del camino
    donde se estudia) para una pareja [{"masa_kg", "o2_L", "kit_kg"?}, {...}].

    Devuelve dict con:
      personas: por astronauta, o2_usado_L, o2_final_pct (lo que queda del tanque al volver), punto_no_retorno
      tanque_pct: por punto, el % de tanque que le quedaría al MÁS justo si volviera desde ahí
      punto_no_retorno: primer punto donde a alguno ya no le alcanza para volver guardando la reserva (o None)
      rescate: margen_L por punto (el peor papel), herido (0/1: quién es el lesionado en el peor caso),
               horas_regreso, sin_rescate (primer punto sin rescate posible, o None) y metros recorridos hasta ahí
      horas: tiempo transcurrido en cada punto
    regresos_calculados: (normal, asistido) de regresos_pareja, para reutilizarlos entre consultas."""
    if len(pareja) != 2:
        raise ValueError("La pareja debe tener exactamente dos astronautas.")
    a = np.asarray(camino)
    idx = (a[:, 0], a[:, 1])
    if radio_m is None:
        radio_m = v3.RES * float(np.hypot(a[:, 0] - base[0], a[:, 1] - base[1]).max()) + 1500.0
    normal, asist = regresos_calculados or regresos_pareja(mapa, base, pareja, radio_m)
    agentes = [astronauta(p.get("masa_kg", 80.0), p.get("kit_kg")) for p in pareja]
    tanques_wh = [p["o2_L"] / v3.O2_L_POR_WH for p in pareja]

    personas, gastado, tanque_pct = [], [], None
    for k, (ag, p) in enumerate(zip(agentes, pareja)):
        pr = v3.perfil_regreso(mapa, camino, {"wh": normal[k]}, ag, tanques_wh[k], reserva, paradas)
        g = np.asarray(pr["gastado_wh"], float)
        gastado.append(g)
        margen_L = np.nan_to_num(np.asarray(pr["margen_wh"], float), nan=-np.inf) * v3.O2_L_POR_WH
        pct = 100 * (margen_L + reserva * p["o2_L"]) / p["o2_L"]      # tanque que quedaría al volver desde ahí
        tanque_pct = pct if tanque_pct is None else np.fmin(tanque_pct, pct)
        personas.append({"masa_kg": float(p.get("masa_kg", 80.0)), "o2_L": float(p["o2_L"]),
                         "o2_usado_L": round(float(g[-1] * v3.O2_L_POR_WH), 1),
                         "o2_final_pct": round(float(100 * (1 - g[-1] * v3.O2_L_POR_WH / p["o2_L"])), 1),
                         "punto_no_retorno": pr["punto_no_retorno"]})

    t_h = tiempos_camino(mapa, camino, agentes[0], paradas)
    h_vuelta = asist["h"][idx]
    margen = np.full(len(a), np.inf)
    herido_peor = np.zeros(len(a), int)
    for herido in (0, 1):
        ayuda = 1 - herido
        mg_ayuda = pareja[ayuda]["o2_L"] - (gastado[ayuda] + asist[f"ayuda{ayuda}"][idx]) * v3.O2_L_POR_WH
        mg_herido = pareja[herido]["o2_L"] - (gastado[herido] + asist[f"herido{herido}"][idx]) * v3.O2_L_POR_WH
        mg = np.fmin(mg_ayuda, mg_herido)
        peor = mg < margen
        margen = np.where(peor, mg, margen)
        herido_peor = np.where(peor, herido, herido_peor)
    margen = np.where(np.isfinite(h_vuelta) & np.isfinite(margen), margen, -1.0)   # sin camino de vuelta
    margen = np.where(t_h + np.nan_to_num(h_vuelta, nan=np.inf) > horas_max, np.minimum(margen, -0.1), margen)
    fuera = np.nonzero(margen < 0)[0]
    sin_rescate = int(fuera[0]) if len(fuera) else None
    no_ret = [p["punto_no_retorno"] for p in personas if p["punto_no_retorno"] is not None]
    d = distancias_camino(mapa, camino)
    return {
        "personas": personas,
        "tanque_pct": np.round(np.nan_to_num(tanque_pct, neginf=-100.0), 1).tolist(),
        "punto_no_retorno": min(no_ret) if no_ret else None,
        "rescate": {
            "margen_L": np.round(margen, 1).tolist(),
            "herido": herido_peor.tolist(),
            "horas_regreso": np.round(np.nan_to_num(h_vuelta, nan=-1.0), 2).tolist(),
            "sin_rescate": sin_rescate,
            "metros_hasta_sin_rescate": None if sin_rescate is None else round(float(d[sin_rescate])),
            "margen_min_L": round(float(margen.min()), 1),
            "supuestos": {"vel_asistida": VEL_ASISTIDA, "carga_ayuda_kg": CARGA_AYUDA_KG, "horas_traje": horas_max},
        },
        "horas": np.round(t_h, 3).tolist(),
        "distancia_m": round(float(d[-1])),
    }


# ---------------------------------------------------------------- radio
def cobertura_radio(mapa, base, radio_m, h_base=ANTENA_HABITAT_M, h_traje=ANTENA_TRAJE_M, paso=2, lote=4000):
    """int8 del tamaño del mapa: 1 = el traje ve la antena del hábitat (línea de vista sobre el relieve, con la
    curvatura de Marte), 0 = sin señal, -1 = no evaluado (fuera de radio_m). Se evalúa una celda de cada `paso`
    (paso 2 = 40 m) y se rellena el resto con la vecina evaluada más cercana."""
    dem = np.nan_to_num(mapa.dem, nan=float(np.nanmin(mapa.dem)))
    n, m = dem.shape
    f0, c0 = base
    k = int(radio_m / v3.RES) + 1
    ff, cc = np.mgrid[max(0, f0 - k):min(n, f0 + k + 1):paso, max(0, c0 - k):min(m, c0 + k + 1):paso]
    forma = ff.shape
    ff, cc = ff.ravel(), cc.ravel()
    D = v3.RES * np.hypot(ff - f0, cc - c0)
    dentro = D <= radio_m
    z0 = dem[f0, c0] + h_base
    ns = max(int(k / paso), 4) + 2                         # ~1 muestra cada `paso` celdas en la línea más larga
    t = np.linspace(0, 1, ns)[1:-1]
    visible = np.zeros(len(ff), bool)
    sel = np.nonzero(dentro)[0]
    for i in range(0, len(sel), lote):
        j = sel[i:i + lote]
        f, c, d = ff[j, None], cc[j, None], D[j, None]
        fs = np.rint(f0 + (f - f0) * t).astype(int)
        cs = np.rint(c0 + (c - c0) * t).astype(int)
        s = d * t
        suelo = dem[fs, cs] + s * (d - s) / (2 * R_MARTE)          # la curvatura levanta el suelo intermedio
        z1 = dem[f[:, 0], c[:, 0]][:, None] + h_traje
        visible[j] = np.all(suelo <= z0 + (z1 - z0) * t, axis=1)
    rejilla = np.where(dentro, visible.astype(np.int8), np.int8(-1)).reshape(forma)
    out = np.full((n, m), -1, np.int8)
    f_lo, c_lo = max(0, f0 - k), max(0, c0 - k)
    bloque = np.repeat(np.repeat(rejilla, paso, axis=0), paso, axis=1)
    alto, ancho = min(n, f0 + k + 1) - f_lo, min(m, c0 + k + 1) - c_lo
    out[f_lo:f_lo + alto, c_lo:c_lo + ancho] = bloque[:alto, :ancho]
    out[f0, c0] = 1
    return out


def con_radio(mapa, cobertura, w=W_SIN_SENAL):
    """Copia del mapa en la que las celdas sin señal (cobertura == 0) cuestan w más. No toca lo no evaluado."""
    return replace(mapa, costo=np.where(cobertura == 0, mapa.costo + w, mapa.costo))


def senal_en_camino(camino, cobertura, mapa=None):
    """Por punto: hay señal (True), no (False) o no se evaluó (None). Y el % del recorrido sin señal y el tramo
    seguido más largo sin señal [m] (lo no evaluado no cuenta como sin señal)."""
    a = np.asarray(camino)
    c = cobertura[a[:, 0], a[:, 1]]
    paso = (np.diff(distancias_camino(mapa, camino)) if mapa is not None
            else v3.RES * np.hypot(*np.diff(a, axis=0).T))
    paso = np.r_[0.0, paso]
    largo = actual = 0.0
    for s, p in zip(c, paso):
        actual = actual + p if s == 0 else 0.0
        largo = max(largo, actual)
    total = float(paso.sum())
    return {"senal": [None if s < 0 else bool(s) for s in c],
            "sin_senal_m": round(float(paso[c == 0].sum())),
            "sin_senal_pct": round(100 * float(paso[c == 0].sum()) / max(total, 1e-9), 1),
            "sin_senal_max_m": round(largo),
            "antenas_m": {"habitat": ANTENA_HABITAT_M, "traje": ANTENA_TRAJE_M}}
