import numpy as np
from quspin.operators import hamiltonian, exp_op
from quspin.basis import spin_basis_1d
from quspin.tools.measurements import obs_vs_time

# Define the Hamiltonian
L = 4 # number of sites
J = 1.0 # interaction strength

# define the spin-1/2 basis
basis = spin_basis_1d(L,pauli=False)

# define the interaction term
J_zz = [[J,i,i+1] for i in range(L-1)] # interaction only between neighbouring sites

# define the Hamiltonian
H = hamiltonian(J_zz, [], basis=basis)

# define the time evolution operator
U = exp_op(-1j*H,'dense')

# define the initial state
init_state = np.random.rand(basis.Ns)

# perform the SSE simulation
times, obs_list = obs_vs_time(U, init_state, return_state=False)

# calculate the expectation value of the energy
energy = obs_list[0]

# plot the energy
import matplotlib.pyplot as plt
plt.plot(times, energy)
plt.xlabel('time')
plt.ylabel('energy')
plt.show()
