import numpy as np
import matplotlib.pyplot as plt

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
VFB = -0.7 # Flatband Voltage
VCB_values = np.arange(0, 1.26, 0.25) #Back-Bias Voltage (VCB)

# Calculation of Fermi Potential, Cox, and Body Factor
phi_f = phi_t * np.log(NA / ni)
Cox = epsilon_ox / tox
gamma = np.sqrt(2 * q * epsilon_s * NA) / Cox

# Define psi_s range
psi_s_values = []
Qc_values = []
VGB_values = []
# Calculate VGB using the given formula
for i in range(len(VCB_values)):
    psi_s_values.append(np.linspace(-0.1,1 + VCB_values[i], 1000))
    Qc = - np.sign(psi_s_values[i]) * np.sqrt(2*q*epsilon_s*NA) * np.sqrt(phi_t*np.exp(-psi_s_values[i]/phi_t) + psi_s_values[i] - phi_t + 
                                                            np.exp(-(2*phi_f+VCB_values[i])/phi_t)*(phi_t*np.exp(psi_s_values[i]/phi_t) -
                                                                                          psi_s_values[i] - phi_t))
    
    VGB = VFB + psi_s_values[i] + np.sign(psi_s_values[i]) * gamma * np.sqrt(
        phi_t * np.exp(-psi_s_values[i] / phi_t) + psi_s_values[i] - phi_t +
        np.exp(-(2 * phi_f+VCB_values[i]) / phi_t) * (phi_t * np.exp(psi_s_values[i] / phi_t) - psi_s_values[i] - phi_t))
    

    Qc_values.append(Qc)
    VGB_values.append(VGB)
    

# Plot VGB vs psi_s for different VCB values
plt.figure(figsize=(10, 6))
for i in range(len(VCB_values)):
    plt.plot(VGB_values[i], Qc_values[i], label=f'VCB = {VCB_values[i]} V')
plt.xlabel('Gate Voltage (VGB, V)')
plt.ylabel('Charge Density (Qs, C/cm^2)')
plt.title('Gate Voltage vs. Charge Density for Different Back-Bias Voltage (VCB)')
plt.legend()
plt.grid(True)
plt.show()

