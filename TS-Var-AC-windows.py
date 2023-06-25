# 2.

# importo paquetes que voy a usar
from cmath import sqrt
import numpy as np
import scipy.stats as st
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Parameters
k = 5 * 10 ** (-7)  # 1/s (Hydraulic constant)
alpha = 1.7 * 10 ** (-4)  # 1/Kelvin (Thermal expansion coefficient)
beta = 8.0 * 10 ** (-4)  # 1/psu (Haline expansion coefficient)
gamma = 2.0 * 10 ** (-8)  # 1/s (temperature relaxation timescale)
V = 10 ** 17  # V = 10**(11)   m^3 (Typical volume of boxes, in a million cubic meters)
S_0 = 35  # psu (reference salinity)
T_star1 = 10  # + 273.15  # Kelvin (Steady sate Tª of box 1)
T_star2 = 15  # + 273.15  # Kelvin


# SOLUTION WITHOUT NOISE
def stommel_trans(T10, T20, S10, S20, T_star1, T_star2, N_steps):
    'Programo cómo resolver la ecuación de Stommel (time-dependent solutions) SIN RUIDO'

    dt = 0.5 * 365 * 24 * 60 * 60  # 1 año en segundos por 0.5

    # Establezco condiciones iniciales para las temperaturas y la salinidad
    T1 = T10  # I set the initial conditions as lists with a single element
    T2 = T20  # Da igual si es Kelvin o ºC (todas in same units) porque sólo calculamos las diferencias de Tª
    S1 = S10
    S2 = S20

    # Creo 'espacios' vacíos del mismo tamaño que el loop (N_steps)
    q = np.empty(N_steps)
    # Varío el Freshwater Flux de manera que tenga el mismo tamaño que el loop y le añado ruido
    F = np.linspace(-0.5, 0.5, N_steps)
    FWF = np.empty(N_steps)  # Freshwater Flux

    for i in range(N_steps):
        # Ecuaciones interdependientes, se van actualizando sus valores a medida que itero en el loop
        q[i] = k * (alpha * (T2 - T1) - beta * (S2 - S1))
        T1 += (np.abs(q[i]) * (T2 - T1) + gamma * (T_star1 - T1)) * dt
        T2 += (np.abs(q[i]) * (T1 - T2) + gamma * (T_star2 - T2)) * dt
        # H = F * S_0 / V, en 1/(psu*s), por lo que F tiene que estar en (m^3/s)
        # De ahí que al definir las constantes se multiplique F*10**6 (ya que 1Sv = 10**6 m^3/s)
        S1 += (np.abs(q[i]) * (S2 - S1) - F[i] * 10 ** 6 * S_0 / V) * dt
        S2 += (np.abs(q[i]) * (S1 - S2) + F[i] * 10 ** 6 * S_0 / V) * dt

        FWF[i] = F[i]  # Almaceno los datos del Freshwater flux

    time = np.arange(N_steps) * dt  # Creo un array para el tiempo del mismo tamaño que el loop

    # Quiero que la función me devuelva el valor de m en Sv (1Sv = 10**6 m^3/s)
    # Donde, m = q (1/s) * V (m^3) / 10**6; he cambiado de unidades
    return q * V / 10 ** 6, FWF, time


# Initial conditions for a strong MOC
T10 = 9.5  # ºC
T20 = 15
S10 = 35  # psu; For a STRONG AMOC, the 'S20-S10' has to be LESS NEGATIVE
S20 = 35.1  #
N_steps = 10 ** 5

m_deterministic, FWF1, time = stommel_trans(T10, T20, S10, S20, T_star1, T_star2, N_steps)
# Almaceno las soluciones del modelo para poder acceder a ellas
time_year = time / 365 / 24 / 3600  # Paso a años la solución


# SOLUTION WITH NOISE
def stommel_trans_noise(T10, T20, S10, S20, T_star1, T_star2, N_steps):
    'Programo cómo resolver la ecuación de Stommel (time-dependent solutions) CON RUIDO'

    dt = 0.5 * 365 * 24 * 60 * 60  # 1 año en segundos por 0.5

    # Establezco condiciones iniciales para las temperaturas y la salinidad
    T1 = T10  # I set the initial conditions as lists with a single element
    T2 = T20  # Da igual si es Kelvin o ºC (todas in same units) porque sólo calculamos las diferencias de Tª
    S1 = S10
    S2 = S20

    # Creo 'espacios' vacíos del mismo tamaño que el loop (N_steps)
    q = np.empty(N_steps)
    seed_value = 13980  # Choose any seed value
    np.random.seed(seed_value)
    # Varío el Freswater Flux de manera que tenga el mismo tamaño que el loop y le añado ruido
    F = np.linspace(-0.5, 0.5, N_steps) + np.random.normal(0, 0.5, N_steps)  # mean 0, standard deviation 0.5
    FWF = np.empty(N_steps)  # Freshwater Flux

    for i in range(N_steps):
        # Ecuaciones interdependientes, se van actualizando sus valores a medida que itero en el loop
        q[i] = k * (alpha * (T2 - T1) - beta * (S2 - S1))
        T1 += (np.abs(q[i]) * (T2 - T1) + gamma * (T_star1 - T1)) * dt
        T2 += (np.abs(q[i]) * (T1 - T2) + gamma * (T_star2 - T2)) * dt
        # H = F * S_0 / V, en 1/(psu*s), por lo que F tiene que estar en (m^3/s)
        # De ahí que al definir las constantes se multiplique F*10**6 (ya que 1Sv = 10**6 m^3/s)
        S1 += (np.abs(q[i]) * (S2 - S1) - F[i] * 10 ** 6 * S_0 / V) * dt
        S2 += (np.abs(q[i]) * (S1 - S2) + F[i] * 10 ** 6 * S_0 / V) * dt

        FWF[i] = F[i]  # Almaceno los datos del Freshwater flux

    time = np.arange(N_steps) * dt  # Creo un array para el tiempo del mismo tamaño que el loop

    # Quiero que la función me devuelva el valor de m en Sv (1Sv = 10**6 m^3/s)
    # Donde, m = q (1/s) * V (m^3) / 10**6; he cambiado de unidades
    return q * V / 10 ** 6, FWF, time


# Initial conditions for a strong MOC
T10 = 9.5  # ºC
T20 = 15
S10 = 35  # psu; For a STRONG AMOC, the 'S20-S10' has to be LESS NEGATIVE
S20 = 35.1  #
N_steps = 10 ** 5

m, FWF2, time2 = stommel_trans_noise(T10, T20, S10, S20, T_star1, T_star2, N_steps)
# Almaceno las soluciones del modelo para poder acceder a ellas
time_year2 = time2 / 365 / 24 / 3600  # / 0.5  # Paso a años la solución


# VARIANCE (cogido de BOERS)

def runstd(x, w):
    n = x.shape[0]
    xs = np.zeros_like(x)

    for i in range(w // 2):
        xw = x[: i + w // 2 + 1]
        xw = xw - xw.mean()

        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1
            xs[i] = np.std(xw)

        else:
            xs[i] = np.nan

    for i in range(n - w // 2, n):
        xw = x[i - w // 2 + 1:]
        xw = xw - xw.mean()

        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1
            xs[i] = np.std(xw)

        else:
            xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2: i + w // 2 + 1]
        xw = xw - xw.mean()

        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1
            xs[i] = np.std(xw)

        else:
            xs[i] = np.nan

    return xs


# AUTO-CORRELATION (cogido de BOERS)
def runac(x, w):
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2):
        xw = x[: i + w // 2 + 1]
        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1

            xs[i] = np.corrcoef(xw[1:], xw[:-1])[0, 1]
        else:
            xs[i] = np.nan

    for i in range(n - w // 2, n):
        xw = x[i - w // 2 + 1:]
        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]

            xw = xw - p0 * np.arange(xw.shape[0]) - p1

            xs[i] = np.corrcoef(xw[1:], xw[:-1])[0, 1]
        else:
            xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2: i + w // 2 + 1]
        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]

            xw = xw - p0 * np.arange(xw.shape[0]) - p1
            xs[i] = np.corrcoef(xw[1:], xw[:-1])[0, 1]
        else:
            xs[i] = np.nan

    return xs


m_std = runstd(m[100:], 2000)
m_ac = runac(m[100:], 2000)

# Create a new FIGURE and SUBPLOTS
fig = plt.figure()
# Add a general title to the entire figure
fig.suptitle('Model with critical transition (AMOC)', fontsize=14, y=0.95)
# First subplot
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(time_year2[100:], m[100:], color="black", label="Stochastic TS")
ax1.plot(time_year[100:], m_deterministic[100:], color="red", label="Deterministic TS")
ax1.set_ylabel('AMOC (Sv)', fontsize=12)
ax1.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax1.set_xlim([0, 50000])  # Set x-axis limits
ax1.legend()
ax1.grid(True)  # Enable grid
# Second subplot
ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(time_year[:78000], m_std[:78000] ** 2, color="darkorange")
ax2.set_ylabel('Variance', fontsize=12)
ax2.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax2.set_xlim(0, 50000)  # Set x-axis limits
ax2.grid(True)  # Enable grid
# Third subplot
ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(time_year[:78000], m_ac[:78000], color="forestgreen")
ax3.set_xlabel('Time (years)', fontsize=12)
ax3.set_ylabel('AC', fontsize=12)
ax3.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax3.set_xlim(0, 50000)  # Set x-axis limits
ax3.grid(True)  # Enable grid
# Add labels above the plots
fig.text(0.05, 0.90, "(a)", fontweight='bold')
fig.text(0.05, 0.62, "(b)", fontweight='bold')
fig.text(0.05, 0.33, "(c)", fontweight='bold')
# Adjust spacing between subplots
fig.tight_layout()
plt.show()

# OTHER WINDOWS
m_std1 = runstd(m[100:], 100)
m_std2 = runstd(m[100:], 1000)
m_std3 = runstd(m[100:], 4000)

m_ac1 = runac(m[100:], 100)
m_ac2 = runac(m[100:], 1000)
m_ac3 = runac(m[100:], 4000)

# Create a new FIGURE and SUBPLOTS
fig = plt.figure()
# Add a general title to the entire figure
fig.suptitle('Model with critical transition (AMOC)', fontsize=14, y=0.95)
# First subplot
ax1 = fig.add_subplot(3, 2, 1)
ax1.plot(time_year[:78000], m_std1[:78000] ** 2, color="darkorange")
ax1.set_title("Variance", fontsize=11)
ax1.set_ylabel('ws = 100', fontsize=11)
ax1.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax1.set_xlim(0, 50000)  # Set x-axis limits
ax1.grid(True)  # Enable grid
# Second subplot
ax2 = fig.add_subplot(3, 2, 3)
ax2.plot(time_year[:78000], m_std2[:78000] ** 2, color="darkorange")
ax2.set_ylabel('ws = 1000', fontsize=11)
ax2.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax2.set_xlim(0, 50000)  # Set x-axis limits
ax2.grid(True)  # Enable grid
# Third subplot
ax3 = fig.add_subplot(3, 2, 5)
ax3.plot(time_year[:78000], m_std3[:78000] ** 2, color="darkorange")
ax3.set_xlabel('Time (years)', fontsize=11)
ax3.set_ylabel('ws = 4000', fontsize=11)
ax3.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax3.set_xlim(0, 50000)  # Set x-axis limits
ax3.grid(True)  # Enable grid

# Fourth subplot
ax4 = fig.add_subplot(3, 2, 2)
ax4.plot(time_year[:78000], m_ac1[:78000], color="forestgreen")
ax4.set_title("Autocorrelation", fontsize=11)
ax4.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax4.set_xlim(0, 50000)  # Set x-axis limits
ax4.grid(True)  # Enable grid
# Fifth subplot
ax5 = fig.add_subplot(3, 2, 4)
ax5.plot(time_year[:78000], m_ac2[:78000], color="forestgreen")
ax5.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax5.set_xlim(0, 50000)  # Set x-axis limits
ax5.grid(True)  # Enable grid
# Sixth subplot
ax6 = fig.add_subplot(3, 2, 6)
ax6.plot(time_year[:78000], m_ac3[:78000], color="forestgreen")
ax6.set_xlabel('Time (years)', fontsize=11)
ax6.axvline(x=40100, color='blue')  # Add x-axis vertical line
ax6.set_xlim(0, 50000)  # Set x-axis limits
ax6.grid(True)  # Enable grid
# Add labels above the plots
fig.text(0.05, 0.90, "(a)", fontweight='bold')
fig.text(0.53, 0.90, "(b)", fontweight='bold')
# Adjust spacing between subplots
fig.tight_layout()
plt.show()
