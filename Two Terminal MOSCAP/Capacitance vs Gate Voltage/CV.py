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

# Calculation of Fermi Potential, Cox, and Body Factor
phi_f = phi_t * np.log(NA / ni)
Cox = epsilon_ox / tox
gamma = np.sqrt(2 * q * epsilon_s * NA) / Cox

# Define psi_s range
psi_s = np.linspace(-0.18, 1.05, 1000)

# Calculate VGB using the given formula
VGB = psi_s + np.sign(psi_s) * gamma *  np.sqrt(
    phi_t * np.exp(-psi_s / phi_t) + psi_s - phi_t +
        np.exp(-2 * phi_f / phi_t) * (phi_t * np.exp(psi_s / phi_t) - psi_s - phi_t))

# Calculate Cgg
Cs = np.where(psi_s == 0, np.sqrt(q * epsilon_s * NA * (1 + np.exp(-2 * phi_f / phi_t)) / phi_t),
              np.sign(psi_s) * np.sqrt(2 * q * epsilon_s * NA) * (
                      1 - np.exp(-psi_s / phi_t) + np.exp(-2 * phi_f / phi_t) * (
                      np.exp(psi_s / phi_t) - 1)) / (2 * np.sqrt(phi_t * np.exp(-psi_s / phi_t) + psi_s - phi_t +
                                                              np.exp(-2 * phi_f / phi_t) * (
                                                                      phi_t * np.exp(psi_s / phi_t) - psi_s - phi_t))))

Cgg = (Cs * Cox)/(Cs + Cox)

# Plot Cgg vs psi_s
plt.figure(figsize=(8, 6))
plt.plot(VGB, Cgg, linewidth=2)
plt.xlabel('Gate Voltage (VGB, V)')
plt.ylabel('MOSCAP Gate Capacitance (Cgg, F/cm^2)')
plt.title('MOSCAP Gate Capacitance vs VGB')
plt.grid(True)
plt.show()
