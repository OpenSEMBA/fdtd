import json
import copy

# Cargar la plantilla
with open('fdtd/test/system/init_solver.fdtd.json', 'r') as f:
    plantilla = json.load(f)

# Defino la malla (50x50x50 celdas de 1cm = Cubo de 50cm)
plantilla["mesh"]["grid"]["numberOfCells"] = [50, 50, 50]
plantilla["mesh"]["grid"]["steps"] = {"x": [0.01], "y": [0.01], "z": [0.01]}

# Defino las coordenadas del Dipolo (en el centro X=25, Y=25)
# El dipolo irá desde Z=5 hasta Z=45 
# El plano de tierra estará en Z=0
plantilla["mesh"]["coordinates"] = [
    {"id": 1, "relativePosition": [25, 25, 5]},  # Base del dipolo (5cm sobre el suelo)
    {"id": 2, "relativePosition": [25, 25, 45]}  # Techo del dipolo
]

# Defino el elemento como una polilínea (el cable)
plantilla["mesh"]["elements"] = [
    {"id": 1, "type": "polyline", "coordinateIds": [1, 2]}
]


# CASO A: FRONTERA MUR + PLANO DE TIERRA
caso_mur = copy.deepcopy(plantilla)
# Configuro todas las paredes como MUR...
caso_mur["boundary"]["all"] = {"type": "mur"}
# la cara inferior (zmin) es el Plano de Tierra (PEC)
caso_mur["boundary"]["zmin"] = {"type": "pec"}

with open('dipolo_mur.json', 'w') as f:
    json.dump(caso_mur, f, indent=4)

# CASO B: FRONTERA PML + PLANO DE TIERRA
caso_pml = copy.deepcopy(plantilla)
# Configuro todas las paredes como PML (absorbentes)...
caso_pml["boundary"]["all"] = {"type": "pml"}
# la cara inferior (zmin) sigue siendo el Plano de Tierra (PEC)
caso_pml["boundary"]["zmin"] = {"type": "pec"}

with open('dipolo_pml.json', 'w') as f:
    json.dump(caso_pml, f, indent=4)

print("Archivos creados con Plano de Tierra (zmin=PEC):")
print("- 'dipolo_mur.json'")
print("- 'dipolo_pml.json'")