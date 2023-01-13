import numpy as np
from scipy.sparse.linalg import eigsh

def Lanczos_ED(H, n_eigvals):
    """
    Perform exact diagonalization of a user-specified Hamiltonian matrix using the Lanczos algorithm.

    Parameters:
    H (ndarray): Hamiltonian matrix
    n_eigvals (int): number of eigenvalues to be computed
    
    Returns:
    eigvals (ndarray): eigenvalues
    eigvecs (ndarray): eigenvectors
    """
    # Obtain the number of states
    n_states = H.shape[0]
    # initialize the Lanczos vectors
    v = np.random.rand(n_states)
    v = v / np.linalg.norm(v)
    w = np.dot(H, v)
    alpha = np.dot(v, w)
    w = w - alpha * v
    beta = np.linalg.norm(w)
    # initialize the tridiagonal matrix
    T = np.zeros((n_eigvals, n_eigvals))
    T[0, 0] = alpha
    # perform Lanczos iterations
    for j in range(1, n_eigvals):
        v = w / beta
        w = np.dot(H, v)
        alpha = np.dot(v, w)
        w = w - alpha * v - beta * v
        beta = np.linalg.norm(w)
        T[j-1, j] = beta
        T[j, j-1] = beta
        T[j, j] = alpha
    # obtain the eigenvalues and eigenvectors of the tridiagonal matrix
    eigvals, eigvecs = np.linalg.eig(T)
    return eigvals, eigvecs

# Example usage:

# Define the Hamiltonian matrix
H = np.array([[1, 2, 3], [2, 4, 5], [3, 5, 6]])

# number of eigenvalues to be computed
n_eigvals = 3

# Perform exact diagonalization
eigvals, eigvecs = Lanczos_ED(H, n_eigvals)

# Print the eigenvalues
print(eigvals)
