import numpy as np
import quspin

# Define the Hamiltonian
L = 2
J = 1.0
t = 0.5

# spin-1/2 operators
Sx, Sy, Sz = quspin.operators.spin(1/2, N=L)

# create spin-1/2 operators for each site
Sx_list = [Sx for i in range(L)]
Sy_list = [Sy for i in range(L)]
Sz_list = [Sz for i in range(L)]

# define the hopping term
hop_x = [[-t, i, (i+1)%L] for i in range(L)]

# define the exchange term
exchange = [[ J, i, (i+1)%L] for i in range(L)]

# create the Hamiltonian
H = quspin.hamiltonians.spin_Hamiltonian(Sx_list, Sy_list, Sz_list, [hop_x, exchange], [], basis=quspin.basis.spinless_fermion_basis_1d)

# diagonalize the Hamiltonian
E, V = H.eigsh(k=4, which='SA')

# print the ground state energy
print(E[0])
