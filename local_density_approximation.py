import numpy as np
from pyscf import gto, dft

# Define the molecule
mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='sto-3g')

# Perform LDA calculation
lda = dft.RKS(mol)
lda.xc = 'LDA'
lda.scf()

# Print the total energy
print(lda.e_tot)

# Obtain the molecular orbitals coefficients and energies
mo_coeff = lda.mo_coeff
mo_energy = lda.mo_energy

# Obtain the electron density
density = lda.make_rdm1()

# Obtain the nuclear-electron attraction energy 
nuclear_repulsion_energy = mol.energy_nuc()

# Obtain the electronic energy
electronic_energy = lda.energy_tot() - nuclear_repulsion_energy
