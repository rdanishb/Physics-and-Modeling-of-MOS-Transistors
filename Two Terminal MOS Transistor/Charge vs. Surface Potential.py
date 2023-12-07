import numpy as np
import matplotlib.pyplot as plt

# Constants
epsilon_0 = 8.854e-14  # Vacuum permittivity in F/cm
q = 1.602e-19  # Elementary charge in C
phi_t = 0.0259  # Thermal voltage in V
ni = 1.5e10 # intrinsic concentration
epsilon_s = 11.7 * epsilon_0 # Silicon permittivity in F/cm

# Given parameters
NA = 1e16  # Acceptor doping concentration in cm^-3
tox = 2e-7  # Oxide thickness in cm (2nm)

# Calculation of Fermi Potential
phi_f = phi_t*(np.log(NA/ni))

# Define surface potential range
psi_s = np.linspace(-0.2, 1, 500) 

# Calculate total semiconductor charge (Qs)
Qs = - np.sign(psi_s) * np.sqrt(2*q*epsilon_s*NA) * np.sqrt(phi_t*np.exp(-psi_s/phi_t) + psi_s - phi_t + 
                                                            np.exp(-2*phi_f/phi_t)*(phi_t*np.exp(psi_s/phi_t) - psi_s - phi_t))

# Calculate inversion charge (Qinv)
Qinv = - np.sqrt(2*q*epsilon_s*NA) * (np.sqrt(psi_s + phi_t*np.exp((psi_s-2*phi_f)/phi_t)) - np.sqrt(psi_s)) 

# Plotting
plt.figure()
plt.plot(psi_s, Qs, 'b', linewidth=2, label='Total Semiconductor Charge')
plt.plot(psi_s, Qinv, 'r', linewidth=2, label='Inversion Charge')
plt.xlabel('Surface Potential (V)')
plt.ylabel('Charge Density (C/cm^2)')
plt.legend()
plt.title('MOSCAP Charge Density vs Surface Potential')
plt.grid(True)
plt.show()
