import numpy as np
from pyscf import gto, scf

# Define the molecule
mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='sto-3g')

# Perform Hartree-Fock calculation
hf = scf.RHF(mol)
hf.scf()

# Print the total energy
print(hf.e_tot)

# Obtain the molecular orbitals coefficients and energies
mo_coeff = hf.mo_coeff
mo_energy = hf.mo_energy

# Obtain the electron density
density = hf.make_rdm1()

# Obtain the nuclear-electron attraction energy 
nuclear_repulsion_energy = mol.energy_nuc()

# Obtain the electronic energy
electronic_energy = hf.energy_tot() - nuclear_repulsion_energy

