import numpy as np

def DMC(num_walkers, num_steps, step_size, V, dt):
    # Initialize the walkers in a random position
    walkers = np.random.uniform(-10, 10, (num_walkers, 3))
    energies = np.zeros(num_walkers)
    for step in range(num_steps):
        # Move the walkers
        walkers += np.random.normal(0, step_size, (num_walkers, 3))
        # Evaluate the energy of each walker
        for i in range(num_walkers):
            energies[i] = V(walkers[i])
        # Calculate the weights of each walker
        weights = np.exp(-dt * energies)
        # Resample the walkers based on their weights
        walkers = np.random.choice(walkers, num_walkers, p=weights/np.sum(weights))
    # Return the final positions of the walkers
    return walkers

# Example usage:

def V(r):
    return r[0]**2 + r[1]**2 + r[2]**2

final_walkers = DMC(1000, 10000, 0.1, V, 0.01)
