import numpy as np
from scipy.optimize import minimize
from scipy.constants import physical_constants
from scipy.integrate import solve_ivp

AU_to_eV = physical_constants["atomic unit of energy"][0] * physical_constants["electron volt"][0]

def dft(r, p, t):
    # Calculate the electronic energy and forces
    E_elec, forces = dft_calculation(r)
    # Calculate the total energy
    E_tot = E_elec + kinetic_energy(p)
    # Calculate the electronic temperature
    T_elec = 2 * E_elec / (3 * len(r))
    # Update the electronic temperature
    t.append(T_elec)
    # Return the derivatives of the positions and momenta
    return p, -forces / (2 * masses)

def kinetic_energy(p):
    return np.sum(p**2 / (2 * masses))

def car_parrinello(r_initial, p_initial, dt, T_elec, T_bath):
    # Optimize the electronic structure
    minimize(dft_calculation, r_initial, args=(), method='BFGS', options={'disp': True})
    # Integrate the equations of motion
    sol = solve_ivp(lambda t, y: dft(y[:len(r_initial)], y[len(r_initial):], t), [0, dt], np.concatenate((r_initial, p_initial)), t_eval=[dt])
    # Update the positions and momenta
    r_initial = sol.y[0][-1]
    p_initial = sol.y[1][-1] + np.random.normal(0, np.sqrt(T_bath * masses / dt), len(r_initial))
    # Scale the momenta to match the electronic temperature
    p_initial *= np.sqrt(T_elec / kinetic_energy(p_initial))
    # Return the updated positions and momenta
    return r_initial, p_initial
