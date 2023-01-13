import numpy as np
import random

def loop_algorithm(L, J, beta, spins):
    """
    Perform the loop algorithm to simulate the Ising model.

    Parameters:
    L (int): number of sites
    J (float): interaction strength
    beta (float): inverse temperature
    spins (ndarray): initial spin configuration

    Returns:
    spins (ndarray): final spin configuration
    """
    # calculate the energy and magnetization
    E, M = energy_and_magnetization(L, J, spins)

    # initialize the loop
    loop = set()
    visited = np.zeros((L, L), dtype=bool)

    # choose a random site
    i = random.randint(0, L - 1)
    j = random.randint(0, L - 1)

    # add the site to the loop
    add_to_loop(i, j, spins, visited, loop, J, beta)

    # perform the loop algorithm
    while len(loop) > 0:
        # choose a random site in the loop
        i, j = random.sample(loop, 1)[0]

        # attempt to add the neighboring sites to the loop
        for ni, nj in neighbors(i, j, L):
            if not visited[ni, nj]:
                add_to_loop(ni, nj, spins, visited, loop, J, beta)

    # update the spin configuration
    for i, j in loop:
        spins[i, j] *= -1

    return spins

def energy_and_magnetization(L, J, spins):
    """
    Calculate the energy and magnetization of the Ising model.

    Parameters:
    L (int): number of sites
    J (float): interaction strength
    spins (ndarray): spin configuration
    Returns:
    E (float): energy
    M (float): magnetization
    """
    E = 0
    M = 0
    for i in range(L):
      for j in range(L):
        E -= J * spins[i, j] * (spins[(i+1)%L, j] + spins[i, (j+1)%L])
        M += spins[i, j]
    return E, M

def add_to_loop(i, j, spins, visited, loop, J, beta):
    """
    Add a site to the loop.
    Parameters:
    i (int): row index of the site
    j (int): column index of the site
    spins (ndarray): spin configuration
    visited (ndarray): visited sites
    loop (set): current loop
    J (float): interaction strength
    beta (float): inverse temperature
    """
    visited[i, j] = True
    loop.add((i, j))
    for ni, nj in neighbors(i, j, L):
        if not visited[ni, nj] and np.exp(-2*beta*J*spins[i,j]*spins[ni,nj]) > random.random():
            add_to_loop(ni, nj, spins, visited, loop, J, beta)
            
def neighbors(i, j, L):
    """
    Get the indices of the neighboring sites.
    Parameters:
    i (int): row index of the site
    j (int): column index of the site
    L (int): number of sites

    Returns:
    neighbors (list): list of neighboring sites
    """
    neighbors = [(i, (j+1)%L), ((i+1)%L, j), (i, (j-1)%L), ((i-1)%L, j)]
    return neighbors

