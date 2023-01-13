import numpy as np

def slater_determinant(wave_functions):
    N = len(wave_functions)
    matrix = np.array([wave_functions[i](x) for i in range(N) for x in range(N)], dtype=np.float64)
    matrix = matrix.reshape((N, N))
    return np.linalg.det(matrix)

def get_wave_functions():
    N = int(input("Enter the number of wave functions: "))
    wave_functions = []
    for i in range(N):
        wave_function = input("Enter a single-particle wave function: ")
        wave_functions.append(lambda x: eval(wave_function))
    return wave_functions

wave_functions = get_wave_functions()
print("The Slater determinant is: ", slater_determinant(wave_functions))
