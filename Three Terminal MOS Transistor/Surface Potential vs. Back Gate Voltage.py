import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
epsilon_0 = 8.854e-14  # Vacuum permittivity in F/cm
q = 1.602e-19  # Elementary charge in C
phi_t = 0.0259  # Thermal voltage in V
ni = 1.5e10  # Intrinsic concentration in cm^-3
epsilon_s = 11.7 * epsilon_0 # Silicon permittivity in F/cm
epsilon_ox = 3.9 * epsilon_0 # SiO2 permittivity in F/cm

# Given parameters
NA = 2e17  # Acceptor doping concentration in cm^-3
tox = 5e-7  # Oxide thickness in cm (5nm)
VFB = -0.7 #Flatband Voltage in V
VGB = 1.5 # Gate Voltage (VCB)


# Calculation of Fermi Potential, Cox, and Body Factor
phi_f = phi_t * np.log(NA / ni)
Cox = epsilon_ox / tox
gamma = np.sqrt(2 * q * epsilon_s * NA) / Cox

#Calculation of VCB
VCB = np.linspace(0, 3, 500)

psi_sa = (-gamma/2 + np.sqrt(((gamma**2)/4) + VGB -VFB))**2

VU = psi_sa - phi_f
VW = psi_sa - 2*phi_f

#VGB = VFB + psi_s + gamma * np.sqrt(psi_s + phi_t * np.exp((psi_s - (2*phi_f + VCB))/phi_t))

# Function to solve
def equation(psi_s, VGB, phi_f, gamma, phi_t, VCB):
    return VGB - (VFB + psi_s + gamma * np.sqrt(psi_s + phi_t * np.exp((psi_s - (2 * phi_f + VCB)) / phi_t)))

# Initial guess for psi_s
initial_guess = 0.1

# Solve the equation for each value of VCB
psi_s_values = []

# Function to solve
def equation(psi_s, VGB, phi_f, gamma, phi_t, VCB):
    return VGB - (VFB + psi_s + gamma * np.sqrt(psi_s + phi_t * np.exp((psi_s - (2 * phi_f + VCB)) / phi_t)))

# Solve the equation for each value of VCB
psi_s_values = []
for vcb_value in VCB:
    if vcb_value < VU:
        # Generate psi from the equation
        initial_guess = 0.1
        solution = fsolve(equation, initial_guess, args=(VGB, phi_f, gamma, phi_t, vcb_value))
        psi_s_values.append(solution[0])
    else:
        # Fix psi_s to psi_sa
        psi_s_values.append(psi_sa)


# Plotting
plt.plot(VCB, psi_s_values, label = "VGFB = 1.5")
plt.xlabel('Back Gate Voltage (V)')
plt.ylabel('Surface Potential ($\psi_s$) (V)')
plt.title('Surface Potential vs. Back Gate Voltage')
plt.legend()
plt.grid(True)
plt.show()
