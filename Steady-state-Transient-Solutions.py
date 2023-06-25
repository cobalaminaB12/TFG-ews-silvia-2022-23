# STEADY-STATE AND TRANSIENT (time-dependent) SOLUTIONS

# importo paquetes que voy a usar
# from cmath import sqrt
import numpy as np
import matplotlib.pyplot as plt

# Model parameters
k = 5*10**(-7)        # 1/s (Hydraulic constant)
alpha = 1.7*10**(-4)  # 1/Kelvin (Thermal expansion coefficient)
beta = 8.0*10**(-4)   # 1/psu (Haline expansion coefficient)
gamma = 2.0*10**(-8)  # 1/s (temperature relaxation timescale)
V = 10**(17)          # V = 10**(11)   m^3 (Typical volume of boxes, in a million cubic meters)
S_0 = 35              # psu (reference salinity)
T_star1 = 10+273.15   # Kelvin (Steady sate Tª of box 1)
T_star2 = 15+273.15   # Kelvin
DT_star = T_star2 - T_star1


# A) STEADY STATE SOLUTION
def stommel_eqn(H):
    q1 = k * alpha * DT_star * (1/2 + np.sqrt(1/4 - beta * H / (k * (alpha * DT_star) ** 2)))
    q2 = k * alpha * DT_star * (1/2 - np.sqrt(1/4 - beta * H / (k * (alpha * DT_star) ** 2)))
    q3 = k * alpha * DT_star * (1/2 + np.sqrt(1/4 + beta * H / (k * (alpha * DT_star) ** 2)))
    q4 = k * alpha * DT_star * (1/2 - np.sqrt(1/4 + beta * H / (k * (alpha * DT_star) ** 2)))
    return q1, q2, q3, q4

# H critical value
H_critical = (k*(alpha*DT_star)**2)/(4*beta)  # 1.1289062499999999e-10
#print("F_critical =", H_critical*V/S_0, "(Sv)")

F_crit = H_critical*V/S_0  # en m^3/s
F_crit = F_crit/10**6  # en Sverdrup
print("The critical value for the freshwater flux is", H_critical, "(m^3/s), ", "i.e., ", F_crit, "(Sv)")

# FOR H > 0
H_values1 = np.linspace(0, 0.5*10**6*S_0/V, 1000)
q_values1 = stommel_eqn(H_values1)
#print(q_values1)
# Plot the results
#plt.plot(H_values*V/S_0/10**6, q_values[0]*V/10**6, label='q1(H)', color="black")                    # q1
#plt.plot(H_values*V/S_0/10**6, q_values[1]*V/10**6, label='q2(H)', color="black", linestyle="--")    # q2
#plt.plot(H_values*V/S_0/10**6, q_values[3]*V/10**6, label='q4(H)', color="red")                      # q4

# FOR H < 0
H_values2 = np.linspace(-0.5*10**6*S_0/V, 0, 1000)
q_values2 = stommel_eqn(H_values2)
# Plot the results
#plt.plot(H_values*V/S_0/10**6, q_values[0]*V/10**6, color="black")  # label='q1(H)', )                    # q1
#plt.xlabel('FWF (Sv)')
#plt.ylabel('AMOC (Sv)')
#plt.title('Steady-state solution')
#plt.legend()
#plt.grid(True)
#plt.show()



# TRANSIENT SOLUTION
# Fixed values
dt = 50000000 #* 365 * 24 * 60 * 60  # seconds, así estaría integrando respecto a 1 año
# Necesito un valor muy grande de dt para ver algo en la gráfica
N_steps = 1000

# El primer q del loop
# q_initial = k * (alpha * (T2 - T1) + beta * (S2 - S1))
# print(q_initial)

def stommel_trans(F, T10, T20, S10, S20, T_star1, T_star2):
    'Programo cómo resolver la ecuación de Stommel (time-dependent solutions)'

    # Establezco condiciones iniciales para las temperaturas y la salinidad
    T1 = T10  # Lo reescribo para que en las ecuaciones quede más compacto y por como están escritas en la teoría
    T2 = T20  # Da igual si es Kelvin o ºC (todas in same units) porque sólo calculamos las diferencias de Tª
    S1 = S10
    S2 = S20

    # I create an empty list of 'q' to start the iterations and accumulate them
    q = []

    for i in range(0, N_steps):
        # Ecuaciones interdependientes, se van actualizando sus valores a medida que itero en el loop
        q.append(k * (alpha * (T2[i] - T1[i]) - beta * (S2[i] - S1[i])))
        T1.append((abs(q[i]) * (T2[i] - T1[i]) + gamma * (T_star1 - T1[i])) * dt + T1[i])
        T2.append((abs(q[i]) * (T1[i] - T2[i]) + gamma * (T_star2 - T2[i])) * dt + T2[i])
        # H = F * S_0 / V, en 1/(psu*s), por lo que F tiene que estar en (m^3/s)
        # De ahí que al definir las constantes se multiplique F*10**6 (ya que 1Sv = 10**6 m^3/s)
        S1.append((abs(q[i]) * (S2[i] - S1[i]) - F * 10 ** 6 * S_0 / V) * dt + S1[i])
        S2.append((abs(q[i]) * (S1[i] - S2[i]) + F * 10 ** 6 * S_0 / V) * dt + S2[i])

    # Quiero que la función me devuelva el valor de m en Sv (1Sv = 10**6 m^3/s)
    # Donde, m = q (1/s) * V (m^3) / 10**6; he cambiado de unidades

    return q[-1] * V / 10 ** 6, T1[-1], T2[-1], S1[-1], S2[-1]
    # Me quedo con el último valor del loop, será solución al modelo de Stommel dependiente del tiempo
    # Para ello selecciono la posición [-1]. Es como escribir list[len(list)-1];
    # Python empieza a contar desde 0, si le pongo como input la longitud de la lista me dará error,
    # porque no existe ese valor, se sale fuera del rango


# Basándonos en las soluciones del steady-state, tenemos como condición inicial m0 = 50 Sv;
# (Case: Strong positive MOC, m > 0; Gradually increase F)

# I create a 'for' loop to gradually increase F from -0.5 to 0.5 Sv
# For that purpose, firstly I create empty lists, to store the set of values
FWF1 = []  # Fresh-Water Flux
list1 = []  # gives the solutions m, T1, T2, S1, S2 of the stommel_transient function
m_list1 = []  # MOC

for f in np.linspace(-0.5, 0.5, 500):
    # I can write f in Sv, because my stommel_trans function changes units from Sv to m^3/s
    # Initial conditions for a strong MOC
    T10 = [10]  # Estoy creado el primer elemento de la lista
    T20 = [15]  # En ºC (son diferencias de Tª, por lo que el factor 273.15K se anula al restar)
    S10 = [35]
    S20 = [35.1]

    # I start adding values to the lists of m and FWF for different values of f (freshwater forcing)
    FWF1.append(f)
    #list1.append(stommel_trans(f, T10, T20, S10, S20, T_star1, T_star2))
    m_list1.append(stommel_trans(f, T10, T20, S10, S20, T_star1, T_star2)[0])

# Plot the results
#plt.plot(FWF1, m_list1, color='black', label="Strong AMOC")
#plt.xlabel('FWF (Sv)')
#plt.ylabel('m (Sv)')
#plt.title('Transient solution: Strong positive MOC')
#plt.grid(True)
#plt.show()


# Basándonos en las soluciones del steady-state, tenemos como condición inicial m0 = 50 Sv;
# (Case 2: Negative MOC, m < 0; Gradually decrease F)
# Initial condition: m0 = -50 (negative)
FWF2 = []  # Fresh-Water Flux
list2 = []  # gives the solutions m, T1, T2, S1, S2 of the stommel_transient function
m_list2 = []  # MOC

for f in np.linspace(0.5, -0.5, 500):
    # I can write f in Sv, because my stommel_trans function changes units from Sv to m^3/s
    # Initial conditions for a strong MOC
    T10 = [13]  # Estoy creado el primer elemento de la lista
    T20 = [15]  # En ºC (son diferencias de Tª, por lo que el factor 273.15K se anula al restar)
    S10 = [31]
    S20 = [36]
    # I start adding values to the lists of m and FWF for different values of f (freshwater forcing)
    FWF2.append(f)
    list2.append(stommel_trans(f, T10, T20, S10, S20, T_star1, T_star2))
    m_list2.append(stommel_trans(f, T10, T20, S10, S20, T_star1, T_star2)[0])

# Plot the results
#plt.plot(FWF2, m_list2, color='blue', label="Weak AMOC")
#plt.xlabel('FWF (Sv)')
#plt.ylabel('AMOC (Sv)')
#plt.title('Transient solution without noise')  #: Negative MOC')
#plt.grid(True)
#plt.legend()
#plt.show()

# Create a new FIGURE and SUBPLOTS
fig = plt.figure()
# Add a general title to the entire figure
fig.suptitle('Model with critical transition (AMOC)', fontsize=16)
# First subplot
ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(H_values1*V/S_0/10**6, q_values1[0]*V/10**6, label='q1(H)', color="black")                    # q1
ax1.plot(H_values1*V/S_0/10**6, q_values1[1]*V/10**6, label='q2(H)', color="black", linestyle="--")    # q2
ax1.plot(H_values1*V/S_0/10**6, q_values1[3]*V/10**6, label='q4(H)', color="red")                      # q4
ax1.plot(H_values2*V/S_0/10**6, q_values2[0]*V/10**6, color="black")
ax1.legend(fontsize=12)
ax1.set_title('a) Steady-state solution', fontsize=14)
ax1.set_xlabel('F (Sv)', fontsize=13)
ax1.set_ylabel('AMOC (Sv)', fontsize=13)
ax1.grid(True)  # Enable grid
# Second subplot
ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(FWF1, m_list1, color='black', label="Strong AMOC")
ax2.plot(FWF2[:253], m_list2[:253], color='blue', label="Weak AMOC")
ax2.set_title('b) Transient Solution without noise', fontsize=14)
ax2.set_xlabel('F (Sv)', fontsize=13)
#ax2.set_ylabel('AMOC (Sv)')
ax2.grid(True)
ax2.legend(fontsize=12)
# Adjust spacing between subplots
fig.tight_layout()
plt.show()

