from qiskit import QuantumCircuit, QuantumRegister
from .hrs_mod_add import cc_carry_add_mod_hrs


def c_mult_hrs(n, a, N):
    """Controlled in-place multiplication using the carry adder (with a constant):
    |x,0...0,0> -> |ax mod N,0...0,0>.

    Register sizes:
    1 control
    n top register |x> -> |ax mod N>
    n+1 ancilla |0> -> |0>

    Approach as described in:
    Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication.
    and
    Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.
    """
    control = QuantumRegister(1, "c")
    topreg = QuantumRegister(n, "x")
    botreg = QuantumRegister(n, "res")
    ancreg = QuantumRegister(1, "anc")

    mult = QuantumCircuit(control, topreg, botreg, ancreg, name="CMULT(%d)MOD(%d)" % (a, N))

    mult.append(c_mult_hrs_partial(n, a, N), mult.qubits)  # = |x, ax mod N, 0> (controlled)

    # SWAP
    for i in range(0, n):
        mult.cswap(control[0], topreg[i], botreg[i])

    ainv = pow(a, -1, N)
    mult.append(c_mult_hrs_partial(n, N - ainv, N), mult.qubits)  # = multiplication x - (a^-1*a)x mod N

    return mult


def c_mult_hrs_partial(n, a, N):
    """Controlled out-of-place multiplication using the carry-adder (with a constant):
    |x,0...0,0> -> |x, ax mod N, 0>.

    Register sizes:
    1 control
    n top register |x> -> |x>
    n bottom register (result) |0> -> |ax mod N>
    1 ancilla |0> -> |0>

    Approach as described in:
    Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication.
    and
    Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.
    """

    control = QuantumRegister(1, "c")
    topreg = QuantumRegister(n, "x")  # |x> - |x>
    botreg = QuantumRegister(n, "res")  # |0> -> |ax mod N>
    ancreg = QuantumRegister(1, "anc")  # |0> -> |0>

    mult = QuantumCircuit(control, topreg, botreg, ancreg, name="CMULTO(%d)MOD(%d)" % (a, N))

    # For each bit add a controlled addition (controlled by topreg and global control)
    for i in range(0, n):
        # need n-1 borrowed qubits. At this point only need one bit from the topreg as control
        # so the n-1 other bits can be borrowed
        borrowed_qubits = [qubit for qubit in topreg if qubit != topreg[i]]
        mult.append(cc_carry_add_mod_hrs(n, (a) * (2 ** i) % N, N),
                    [control[0], topreg[i]] + list(botreg) + borrowed_qubits + list(ancreg))

    return mult