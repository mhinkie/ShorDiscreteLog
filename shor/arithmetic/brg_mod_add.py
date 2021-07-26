from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import QFT
from .draper_add import phi_add_builtin, phi_addition_controlled, phi_addition_cc


def phi_add_mod_N_cc(n, y, N):
    """Modular addition in fourier space as described in:
    Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.

    Uses two control qubits.

    |phi(a)>|0> -> |phi((a+y) mod N)>|0>.

    Register sizes:
    2 control
    n |phi(a)> -> |phi((a+y) mod N)>
    1 ancilla
    """
    ctrl = QuantumRegister(2, "ctrl")
    reg = QuantumRegister(n, "a")  # Holds the input
    aux = QuantumRegister(1, "aux")  # Holds the msb after the first step to use as control later

    qc_add = QuantumCircuit(ctrl, reg, aux, name="CCADD(%d)MOD(%d)" % (y, N))

    iqft = QFT(n, do_swaps=False).inverse()
    qft = QFT(n, do_swaps=False)

    # add y
    qc_add.append(phi_addition_cc(n, y, 1), list(ctrl) + list(reg))

    # decrease by N
    qc_add.append(phi_add_builtin(n, N, -1), reg)

    # iqft
    qc_add.append(iqft, reg)

    # msb is in highest register
    qc_add.cx(reg[n - 1], aux[0])

    # QFT again, because we add in fourier space again
    qc_add.append(qft, reg)

    # controlled add
    cadd = phi_addition_controlled(n, N, 1)
    qc_add.append(cadd, [aux[0]] + list(reg))  # have to reorder here so the last bit controls the gate

    # substract y again
    qc_add.append(phi_addition_cc(n, y, -1), list(ctrl) + list(reg))

    # do inv. qft to get msb again
    qc_add.append(iqft, reg)

    # flip msb
    qc_add.x(reg[n - 1])

    # flip aux according to msb
    qc_add.cx(reg[n - 1], aux[0])

    # flip it back
    qc_add.x(reg[n - 1])

    # back to fourier space
    qc_add.append(qft, reg)

    # Now add y back again to get to the proper result again
    qc_add.append(phi_addition_cc(n, y, 1), list(ctrl) + list(reg))

    return qc_add


def phi_add_mod_N(n, y, N):
    """Modular addition in fourier space as described in:
    Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.

    |phi(a)>|0> -> |phi((a+y) mod N)>|0>.

    Register sizes:
    n |phi(a)> -> |phi((a+y) mod N)>
    1 ancilla
    """
    reg = QuantumRegister(n, "a")  # Holds the input
    aux = QuantumRegister(1, "aux")  # Holds the msb after the first step to use as control later
    aux_index = n

    qc_add = QuantumCircuit(reg, aux, name="ADD(%d)MOD(%d)" % (y, N))

    iqft = QFT(n, do_swaps=False).inverse()
    qft = QFT(n, do_swaps=False)

    # add y
    qc_add.append(phi_add_builtin(n, y, 1), reg)

    # decrease by N
    qc_add.append(phi_add_builtin(n, N, -1), reg)

    # iqft
    qc_add.append(iqft, reg)

    # msb is in highest register
    qc_add.cx(reg[ n -1], aux[0])

    # QFT again, because we add in fourier space again
    qc_add.append(qft, reg)

    # controlled add
    cadd = phi_addition_controlled(n, N, 1)
    qc_add.append(cadd, [aux_index]+list(reg))  # have to reorder here so the last bit controls the gate

    # substract y again
    qc_add.append(phi_add_builtin(n, y, -1), reg)

    # do inv. qft to get msb again
    qc_add.append(iqft, reg)

    # flip msb
    qc_add.x(reg[n-1])

    # flip aux according to msb
    qc_add.cx(reg[n-1], aux[0])

    # flip it back
    qc_add.x(reg[n-1])

    # back to fourier space
    qc_add.append(qft, reg)

    # Now add y back again to get to the proper result again
    qc_add.append(phi_add_builtin(n, y, 1), reg)

    return qc_add