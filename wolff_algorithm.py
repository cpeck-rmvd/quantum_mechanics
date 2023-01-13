import numpy as np
import random

def Wolff_algorithm(L, J, beta, spins):
    """
    Perform the Wolff algorithm to simulate the Ising model.

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

    # initialize the cluster
    cluster = set()
    visited = np.zeros((L, L), dtype=bool)

    # choose a random site
    i = random.randint(0, L - 1)
    j = random.randint(0, L - 1)

    # add the site to the cluster
    add_to_cluster(i, j, spins, visited, cluster, J, beta)

    # perform the Wolff algorithm
    while len(cluster) > 0:
        # choose a random site in the cluster
        i, j = random.sample(cluster, 1)[0]

        # attempt to add the neighboring sites to the cluster
        for ni, nj in neighbors(i, j, L):
            if not visited[ni, nj]:
                add_to_cluster(ni, nj, spins, visited, cluster, J, beta)
    # update the spin configuration
    for i, j in cluster:
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

def add_to_cluster(i, j, spins, visited, cluster, J, beta):
"""
Add a site to the cluster.
Parameters:
i (int): row index of the site
j (int): column index of the site
spins (ndarray): spin configuration
visited (ndarray): visited sites
cluster (set): current cluster
J (float): interaction strength
beta (float): inverse temperature
"""
  visited[i, j] = True
  cluster.add((i, j))
  for ni, nj in neighbors(i, j, L):
      if not visited[ni, nj] and np.exp(-2*beta*J*spins[i,j]*spins[ni,nj]) > random.random():
          add_to_cluster(ni, nj, spins, visited, cluster, J, beta)
          
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



