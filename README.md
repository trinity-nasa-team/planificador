# planificador_v3

Planificador de misiones NSGA-II (Jezero) listo para correr por sí solo.

- `mision_v3.py` – el planificador (mapa de costos, rutas normal/segura, NSGA-II, búsqueda local, ILS)
- `exportar_mision_v3.py` – genera `salida/rutas_v3.json` para el visor 3D
- `proto1_capas.py` – carga las capas desde `data/processed/`
- `demo_planificador_v3.ipynb` – notebook paso a paso
- `salida/rutas_v3.json` – ejemplo ya generado (factor 2.0, semilla 42)

## Uso
```
pip install -r requirements.txt     # Python 3.10+ (probado en 3.12)
python exportar_mision_v3.py --factor 2.0 --semilla 42
jupyter notebook demo_planificador_v3.ipynb      # abrir desde esta carpeta
```

## API local (la Raspberry Pi pide las rutas a la laptop)

La laptop corre el planificador como servicio HTTP; la Raspberry solo necesita `requests`.

**En la laptop** (elige una):
```
pip install -r requirements-api.txt
uvicorn api:app --host 0.0.0.0 --port 8000
```
```
docker compose up -d --build
```
Abre http://localhost:8000/docs para probar los endpoints desde el navegador.

Busca la IP de la laptop con `ipconfig` (Dirección IPv4, p. ej. 192.168.1.50) y, la primera vez, permite el puerto en el firewall de Windows (PowerShell como administrador):
```
New-NetFirewallRule -DisplayName "Planificador 8000" -Direction Inbound -Protocol TCP -LocalPort 8000 -Action Allow
```

**En la Raspberry** (misma red Wi-Fi):
```
pip install requests
python cliente_rpi.py --servidor http://192.168.1.50:8000 --prioridad segura
```

| Endpoint | Para qué |
|---|---|
| `GET /salud` | ¿el servidor responde? ¿está ocupado calculando? |
| `GET /ejemplo` | inicio y meta de la ruta de Perseverance para probar |
| `POST /planear` | `{"inicio":[x,y], "meta":[x,y], "prioridad":"segura"}` → waypoints, objetivos, métricas |
| `POST /replanear` | igual, desde la posición actual y con `"visitados":[[x,y],...]` para no repetir objetivos |

Coordenadas en metros locales (x = este, y = norte), las mismas de `rutas_v3.json`. Prioridades: `ciencia`, `rapida`, `segura`, `equilibrada`, `pesos`. Misión de ida y vuelta (EVA): `"meta": null, "horas": 8, "radio_m": 2000, "agente": "astronauta"`.

Un plan tarda ~10 s; uno repetido sale de caché al instante. Si se cae la conexión, el cliente usa `ultimo_plan.json`. Clave opcional: arranca con `API_KEY=algo` y el cliente con `--api-key algo`.
