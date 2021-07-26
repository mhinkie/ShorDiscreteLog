from qiskit import transpile, QuantumCircuit, QuantumRegister
import numpy as np

def set_input(n, value_dec):
    """sets input in lsb first repr"""
    value_bin = np.binary_repr(value_dec, n)
    #start with lsb
    value_bin = value_bin[::-1]

    reg = QuantumRegister(n, "r")
    circ = QuantumCircuit(reg, name="input=%d" % value_dec)

    for i in range(0, n):
        if value_bin[i] == '1':
            circ.x(reg[i])

    return circ.to_instruction()


def run_on_simulator(circ, sim, rep):
    circ = transpile(circ, sim)

    result = sim.run(circ, shots=rep).result()
    counts = result.get_counts(circ)

    return counts
