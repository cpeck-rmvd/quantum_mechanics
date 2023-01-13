import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import execute, Aer

# Define the number of qubits
n_qubits = 3

# Create a quantum circuit
q = QuantumRegister(n_qubits)
c = ClassicalRegister(n_qubits)
circuit = QuantumCircuit(q, c)

# Apply the Hadamard gate to each qubit
for i in range(n_qubits):
    circuit.h(q[i])

# Apply a controlled-Z gate between qubits 0 and 1
circuit.cz(q[0], q[1])

# Measure the qubits
for i in range(n_qubits):
    circuit.measure(q[i], c[i])

# Execute the circuit on a simulator
backend = Aer.get_backend('qasm_simulator')
result = execute(circuit, backend, shots=10000).result()
counts = result.get_counts()

# Print the results
print(counts)
