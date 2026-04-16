import numpy as np
import fdtd
import matplotlib.pyplot as plt

def ejecutar_simulacion(tipo_borde="PML"):
    # Definición de la malla
    grid = fdtd.Grid(shape=(60, 60, 40), grid_spacing=0.01)

    if tipo_borde == "PML":
        grid[0:10, :, :] = fdtd.PML(name="pxlow")
        grid[-10:, :, :] = fdtd.PML(name="pxhigh")
        grid[:, 0:10, :] = fdtd.PML(name="pylow")
        grid[:, -10:, :] = fdtd.PML(name="pyhigh")
        grid[:, :, -10:] = fdtd.PML(name="pzhigh")
    
    # Fuente (Dipolo Z=5 a 7)
    grid[30, 30, 5:7] = fdtd.LineSource(period=1/1e9, name="fuente")

    # Configuración de medida
    R = 15 #radio de medición: 15 celdas de distancia desde la antena
    angulos_deg = np.array([10, 20, 30, 40, 50, 60, 70, 80])
    max_valores = np.zeros(len(angulos_deg))

    print(f"Simulando caso: {tipo_borde}...")
    
    for t in range(800):
        grid.step()
        grid.E[:, :, 0, :] = 0  # Suelo PEC
        
        for i, ang in enumerate(angulos_deg):
            rad = np.radians(ang)
            x = int(30 + R * np.cos(rad))
            z = int(R * np.sin(rad))
            
            # Guardo el valor máximo absoluto detectado
            val_actual = np.abs(float(grid.E[x, 30, z, 2]))
            if val_actual > max_valores[i]:
                max_valores[i] = val_actual

    # Evitar división por cero y normalizar
    if np.max(max_valores) > 0:
        return max_valores / np.max(max_valores)
    return max_valores

# Ejecución

puntos_pml = ejecutar_simulacion("PML")
puntos_mur = ejecutar_simulacion("Mur")

# Teoría

'''
theta_plot = np.linspace(0, 90, 100)
theta_rad = np.radians(theta_plot)
teoria = np.cos(theta_rad)**10  
teoria /= np.max(teoria)
'''

h = 0.06
lambd = 0.3
k = 2 * np.pi / lambd

theta_plot = np.linspace(0, 90, 100)
alpha = np.radians(theta_plot) # elevación desde el suelo

# Patrón del dipolo: cos(alpha) es el equivalente al sin(theta) de Balanis
patron = np.cos(alpha)

# Factor de grupo (Interferencia): cos(kh sin alpha) es el cos(kh cos theta) de Balanis
interferencia = np.abs(np.cos(k * h * np.sin(alpha)))

# Proyección sobre Ez: Como medimos el campo vertical, proyectamos con cos(alpha)
proyeccion_ez = np.cos(alpha)

# Teoría final: Amplitud de Ez
teoria = patron * interferencia * proyeccion_ez
teoria /= np.max(teoria)


plt.figure(figsize=(10, 6))
plt.plot(theta_plot, teoria, 'b-', label='Teoría', linewidth=2)

ang_exp = np.array([10, 20, 30, 40, 50, 60, 70, 80])
plt.scatter(ang_exp, puntos_pml, color='red', s=70, label='FDTD + PML', edgecolors='black', zorder=5)
plt.scatter(ang_exp, puntos_mur, color='green', marker='x', s=100, label='FDTD + Mur', zorder=4)

plt.xlabel("Ángulo de Elevación (grados)")
plt.ylabel("Intensidad Ez Normalizada")
plt.ylim(-0.05, 1.1)
plt.xlim(0, 90)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

plt.savefig("grafica_final.png")
plt.show()
