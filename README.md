# planificador_v3

Planificador de misiones NSGA-II (Jezero) listo para correr por sí solo.

- `mision_v3.py` – el planificador (mapa de costos, rutas normal/segura, NSGA-II, búsqueda local, ILS)
- `exportar_mision_v3.py` – genera `salida/rutas_v3.json` para el visor 3D
- `proto1_capas.py` – carga las capas desde `data/processed/`
- `demo_planificador_v3.ipynb` – notebook paso a paso
- `salida/rutas_v3.json` – ejemplo ya generado (factor 2.0, semilla 42), con coordenadas locales y lon/lat
- `validar_seitah.py` – compara las rutas con el recorrido real de Perseverance (tramos A, B, C)
- `herramientas/preparar_capas_marswalk.py` – regenera las capas de arena y pendiente HiRISE
- `tests/` – 59 pruebas automáticas
- **Cambios de Marswalk v4: ver `CAMBIOS_V4.md`** (arena de Séítah, pendiente HiRISE, IA + robot, astronauta con O₂ y punto de no retorno, correcciones)

## Uso
```
pip install -r requirements.txt     # Python 3.10+ (probado en 3.12)
python exportar_mision_v3.py --factor 2.0 --semilla 42
jupyter notebook demo_planificador_v3.ipynb      # abrir desde esta carpeta
python validar_seitah.py                          # tramo B de Séítah: v3 498 m -> v4 125 m de la ruta real
cd tests && python -m pytest -q                   # 59 pruebas
```
