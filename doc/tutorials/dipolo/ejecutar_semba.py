import json
import fdtd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def generar_animacion(archivo_json):
    print(f"--- Iniciando Animación desde: {archivo_json} ---")
    with open(archivo_json, 'r') as f:
        conf = json.load(f)

    celdas = conf["mesh"]["grid"]["numberOfCells"] 
    paso = conf["mesh"]["grid"]["steps"]["x"][0]   
    
    try:
        grid = fdtd.Grid(shape=(celdas[0], celdas[1], celdas[2]), grid_spacing=paso)
    except:
        grid = fdtd.grid(shape=(celdas[0], celdas[1], celdas[2]), grid_spacing=paso)

    # Fronteras (PML/Mur)
    tipo_borde = conf["boundary"]["all"]["type"].lower()
    if tipo_borde == "pml":
        grid[0:5, :, :] = fdtd.PML(name="pml_x_low")
        grid[-5:, :, :] = fdtd.PML(name="pml_x_high")
        grid[:, 0:5, :] = fdtd.PML(name="pml_y_low")
        grid[:, -5:, :] = fdtd.PML(name="pml_y_high")

    # Dipolo (Holland Model)
    grid[25, 25, 5:45] = fdtd.LineSource(period=15e-10, name="dipolo")

    fig, ax = plt.subplots(figsize=(8, 6))
    fotogramas = []
    
    for t in range(100): # 100 pasos de tiempo
        grid.step()
        # Aplico el Plano de Tierra (PEC) en cada paso
        grid.E[:, :, 0, 2] = 0 
        
        # Capturo un corte vertical (plano XZ) para ver la radiación
        # Accedo a los datos crudos del campo Ez
        data = np.abs(grid.E[:, 25, :, 2]) 
        im = ax.imshow(data.T, cmap='viridis', origin='lower', animated=True)
        if t == 0:
            ax.imshow(data.T, cmap='viridis', origin='lower') # Fondo inicial
        fotogramas.append([im])
        
        if t % 20 == 0:
            print(f"Frame {t}/100...")

    ani = animation.ArtistAnimation(fig, fotogramas, interval=50, blit=True)
    
    nombre_video = archivo_json.replace(".json", ".gif")
    print(f"Guardando animación en {nombre_video}...")
    ani.save(nombre_video, writer='pillow')
    
    plt.title(f"Radiación Dipolo (Corte XZ) - {tipo_borde.upper()}")
    plt.xlabel("X (celdas)")
    plt.ylabel("Z (Altura)")
    plt.show()

if __name__ == "__main__":
    generar_animacion("dipolo_pml.json")
    generar_animacion("dipolo_mur.json")