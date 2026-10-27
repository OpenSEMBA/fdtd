# Experimento SIGSTOP: un rango MPI detenido bloquea a los demás con point probes

Fecha: 2026-10-27
Binario: `build-rls-mpi/bin/semba-fdtd` (commit `a40a22ec`, Release, MPI, 2 rangos)

## Objetivo

Comprobar, sin modificar el código del solver, qué ocurre cuando un proceso MPI se
detiene en mitad de una simulación que solo tiene point probes. El objetivo es
distinguir si el bloqueo que se observa con point probes viene del camino de
salida o de la sincronización del bucle temporal.

## Método

1. Caso: `testData/cases/stability/lossy_cavity.fdtd.json`, copiado a un directorio
   de trabajo. Cambios sobre la copia (no sobre el original):
   - `numberOfSteps` = 10 000 000 (el original tiene 8000 y terminaría en segundos).
   - Un segundo point probe `ez_rank0` en el nodo (4,5,2), que cae en la franja de
     z 0–6 del rango 0. El probe original `ez` está en (13,9,7), en la franja del
     rango 1. Ver `lossy_cavity_two_probes.fdtd.json`.
   - Hace falta `gauss_cavity.exc` junto al caso (está en `testData/cases/stability/`).
2. Lanzamiento: `mpirun -np 2 semba-fdtd -i lossy_cavity.fdtd.json`.
3. Identificación de los PIDs por la variable `OMPI_COMM_WORLD_RANK` del entorno:
   rango 0 = PID 92242, rango 1 = PID 92243 (propietario del probe `ez`).
4. `sigstop_experiment.sh <pid_rango1>` muestrea cada 5 s durante 30 s en tres fases:
   - *baseline*: ambos procesos corriendo.
   - *paused*: `kill -STOP` sobre el rango 1.
   - *resumed*: `kill -CONT` sobre el rango 1.
   En cada muestra registra el estado de cada proceso (`/proc/<pid>/stat`), su tiempo
   de CPU y el número de filas de cada fichero de probe.

## Resultados

Traza completa: `sigstop_timeline.txt`.

| Fase | Rango 1 (`SIGSTOP`) | Filas probe `ez` (rango 1) | Filas probe `ez_rank0` (rango 0) | CPU rango 0 |
|---|---|---|---|---|
| baseline (30 s) | `R` | 757 121 → 1 486 721 | 757 377 → 1 486 977 | `R`, ~100 % |
| paused (30 s) | `T` (parado) | **1 627 393 (sin cambios)** | **1 627 393 (sin cambios)** | `R`, ~100 % |
| resumed (30 s) | `R` | 1 628 033 → 2 392 705 | 1 628 417 → 2 393 089 | `R`, ~100 % |

Observaciones:

- Durante la pausa, **el probe del rango 0 también se congela** en la misma fila que el
  del rango 1. El rango 0 no tiene ninguna relación directa con el rango 1 salvo el
  intercambio de campos de frontera en cada paso.
- El rango 0 sigue al 100 % de CPU durante toda la pausa, sin avanzar. Eso es la espera
  activa de OpenMPI dentro de la comunicación, no progreso. No hay ningún error, mensaje
  ni timeout.
- Antes de la pausa el rango 0 iba unas 256 filas por delante del rango 1. Esa diferencia
  es el desfase de escritura de cada fichero, no un avance real: los dos procesos están
  acoplados en paso.
- Al reanudar, ambos continúan sin intervención y no aparece ningún `FATAL` ni `ABORT`
  en `run.log` (comprobado: 0 coincidencias).

## Interpretación

- **El bloqueo se propaga al siguiente intercambio de frontera.** Con point probes no
  hay ninguna escritura MPI en el camino de salida (`pointProbeOutput.F90` es local).
  El bloqueo viene del bucle temporal: en cada paso hay `MPI_Barrier` y
  `FlushMPI_E_Cray` / `FlushMPI_H_Cray` (`timestepping.F90` ~2049 y ~2086;
  `mpicomm.F90:1425–1593`), con `MPI_WAITALL` sobre los vecinos.
- **Un rango detenido detiene a todos en el siguiente paso.** La cadena de dependencias
  hace que los demás no puedan avanzar más de un paso sin el rango parado.
- **Esto explica el síntoma observado con point probes**, sin depender de movies ni del
  relay de `relay_flush_batch`. Cualquier retraso de un rango (I/O bloqueado, carga
  desigual, un proceso colgado) se traslada a todos los rangos.
- **Sin timeout, un rango parado de forma permanente es un cuelgue permanente.** En
  este experimento la pausa duró 30 s y la simulación se recuperó sola. Si el proceso
  no se reanudara, la simulación quedaría bloqueada indefinidamente sin error.

## Límites del experimento

- Solo se ha pausado el rango 1. No se ha hecho la prueba simétrica (pausar el rango 0),
  aunque por la lockstep de los pasos se espera el mismo resultado.
- Se ha medido el avance con el número de filas de los ficheros de probe, no con un
  contador de pasos. Las filas se escriben en bloques, así que la resolución es de
  pocos miles de pasos.
- No hay `gdb` en el sistema, así que no se ha obtenido la pila de cada proceso.
  La afirmación de que el rango 0 está en `WAITALL` o `Barrier` se deduce del acoplamiento
  y de la espera activa, no de una traza de pila.
- Un proceso `SIGSTOP`-eado no es lo mismo que un proceso que tarda: aquí la detención
  es total y conocida. Un I/O lento o un desequilibrio tendrían una firma distinta
  (avance lento, no detención).

## Cómo reproducirlo

```bash
cp lossy_cavity_two_probes.fdtd.json <dir>/lossy_cavity.fdtd.json
cp testData/cases/stability/gauss_cavity.exc <dir>/
cd <dir>
mpirun -np 2 <repo>/build-rls-mpi/bin/semba-fdtd -i lossy_cavity.fdtd.json > run.log 2>&1 &
# identificar PIDs por OMPI_COMM_WORLD_RANK en /proc/<pid>/environ
bash <repo>/doc/hang-experiments/sigstop_experiment.sh <pid_rango1> > timeline.txt
```

El script asume el nombre de los ficheros de probe de este caso
(`lossy_cavity.fdtd_ez_Ez_13_9_7_tm.dat` y `lossy_cavity.fdtd_ez_rank0_Ez_4_5_2_tm.dat`).
Parar la ejecución al terminar con `pkill -f "[s]emba-fdtd -i"`.
