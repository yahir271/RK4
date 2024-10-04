import numpy as np
import matplotlib.pylab as plt

# Definir el estado para el progreso
state = [1]

def print_percent_done(index, total, state, title='Please wait'):
    percent_done2 = (index + 1) / total * 100
    percent_done = round(percent_done2, 1)

    print(f'\t⏳{title}: {percent_done}% done', end='\r')

    if percent_done2 > 99.9 and state[0]:
        print('\t✅')
        state[0] = 0

# Parámetros iniciales
w0 = 2/3
h = 0.01
t_0 = 0
x_0 = 0.2
y_0 = 0
num_plot_points = 100

# Definir funciones f y g
def f(t, x, y):
    return y #omega

def g(t, x, y):
    q = 6     #periodo4= 0.25,0.75, caotico=1.05 periodo3=0.5,0.358
    return -(y/q) - np.sin(x) + (2) * np.cos((2/3) * t)

# Implementación del paso de Runge-Kutta de cuarto orden
def RK4step(t, x, y, h):
    k1 = f(t, x, y)
    u1 = g(t, x, y)

    k2 = f(t + h/2, x + h*k1/2, y + h*u1/2)
    u2 = g(t + h/2, x + h*k1/2, y + h*u1/2)

    k3 = f(t + h/2, x + h*k2/2, y + h*u2/2)
    u3 = g(t + h/2, x + h*k2/2, y + h*u2/2)

    k4 = f(t + h, x + h*k3, y + h*u3)
    u4 = g(t + h, x + h*k3, y + h*u3)

    x_next = x + (k1 + 2*k2 + 2*k3 + k4) * h / 6
    y_next = y + (u1 + 2*u2 + 2*u3 + u4) * h / 6

    return x_next, y_next

# Inicializar variables para la sección de Poincaré
x_section, y_section = [], []
T = 2 * np.pi / w0  # Periodo del sistema

# Loop para resolver el sistema y graficar la sección de Poincaré
x, y = x_0, y_0
for k in range(num_plot_points):
    t = 0
    while t < T - 2 * h: #1.2
        x, y = RK4step(t, x, y, h)
        t += h
    # Realizar el último paso para llegar exactamente a t = T
    x, y = RK4step(t, x, y, T - t)
    
    if k % 100 == 0:
        print_percent_done(k, num_plot_points, state, 'Solving given DE')
    
    x_section.append(x)
    y_section.append(y)

# Graficar la sección de Poincaré
plt.plot([x % (2 * np.pi) for x in x_section], y_section, 'gs', markersize=2)
plt.title('Sección de Poincaré para q = 6')
plt.xlabel('theta',fontsize =13)
plt.ylabel('$\omega$',fontsize = 13)
plt.show()

