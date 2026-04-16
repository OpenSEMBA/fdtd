import numpy as np
import matplotlib.pyplot as plt

frecuencia = 300e6  # 300 MHz -> lambda = 1m
L = 0.5            # dipolo de media onda (lambda/2)
h = 0.1            # Altura sobre el plano de tierra
c = 3e8
k = 2 * np.pi / (c / frecuencia) # Número de onda

# Ángulos (de 0 a 180 grados, solo la parte de arriba de la tierra)
theta = np.linspace(0.001, np.pi/2, 360) 

# Solución Analítica (Far Field - Campo Lejano) - Referencia 3
# Campo de un dipolo sobre plano de tierra:
term1 = np.abs((np.cos((k*L/2)*np.cos(theta)) - np.cos(k*L/2)) / np.sin(theta))
# Factor de Array debido a la "Antena Imagen" del plano de tierra
term2 = np.abs(np.sin(k * h * np.cos(theta)))

# Patrón total (Normalizado)
patron_total = term1 * term2
patron_total /= np.max(patron_total)

# Gráfica Polar
plt.figure(figsize=(7,7))
ax = plt.subplot(111, projection='polar')
ax.plot(theta, patron_total, label='Solución Analítica', color='red')

# Como es sobre plano de tierra (PEC), solo hay radiación de 0 a 90 grados
# (de 90 a -90 respecto a la normal)
ax.plot(-theta, patron_total, color='red')

ax.set_theta_zero_location('N') # El norte es el eje Z (antena vertical)
ax.set_theta_direction(-1)     # Dirección horaria
ax.set_title("Patrón de Radiación (Far Field)\nDipolo sobre Plano de Tierra")
ax.legend(loc='lower center')
plt.savefig('patron_radiacion_analitico.png')
plt.show()