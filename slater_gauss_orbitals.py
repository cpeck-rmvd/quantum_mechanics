import numpy as np

def sto(n, l, m, r, Z):
    rnl = (2*Z/n)**(3/2) * np.exp(-Z*r/n) * (2*r/n)**l
    return rnl * Ylm(l, m, r)

def gto(n, l, m, r, Z, alpha):
    rnl = (2*alpha/np.pi)**(3/4) * (alpha*r)**l * np.exp(-alpha*r**2/2)
    return rnl * Ylm(l, m, r)

# Ylm spherical harmonics function
def Ylm(l, m, r):
    theta, phi = cartesian_to_spherical(r)
    return special.sph_harm(m, l, phi, theta)

# Example usage:

# Nuclear charge
Z = 1

# STO parameters
n = 2
l = 0
m = 0

# GTO parameters
alpha = 1

# Grid for plotting the orbitals
r = np.linspace(0, 10, 100)

# Calculate the STO and GTO orbitals
sto_values = sto(n, l, m, r, Z)
gto_values = gto(n, l, m, r, Z, alpha)

# Plot the STO and GTO orbitals
import matplotlib.pyplot as plt
plt.plot(r, sto_values, label='STO')
plt.plot(r, gto_values, label='GTO')
plt.legend()
plt.show()
