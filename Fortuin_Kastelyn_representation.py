import networkx as nx
from networkx.algorithms.approximation import independent_set
import random

# Define the Ising model graph
L = 4 # number of sites
J = 1.0 # interaction strength
beta = 1.0 # inverse temperature

# create a 2D lattice graph
G = nx.grid_2d_graph(L, L, periodic=True)

# assign weights to the edges
for (u, v) in G.edges():
    G[u][v]['weight'] = -J * (random.random() < 0.5)

# find the independent sets
ind_sets = independent_set(G, weight='weight')

# calculate the partition function
Z = sum(np.exp(-beta*sum(G[u][v]['weight'] for (u, v) in ind_set)) for ind_set in ind_sets)

