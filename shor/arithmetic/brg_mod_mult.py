from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import QFT
from .brg_mod_add import phi_add_mod_N_cc
import math


def mult_mod_N_c(a, N, size_requested=None):
    """Controlled in-place multiplication:
    |x>|0> -> |ax mod N>|0>

    As presented in: Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.

    Register sizes:
    1 control
    n input=output |x> -> |ax mod N>
    n+2 ancilla

    n can either be supplied (size_requested) or will be calulated as ceil(ld(N)).
    The module N should be co-prime to the input constant a."""
    size = 0
    if size_requested == None:
        size = math.ceil(math.log(N, 2))
    else:
        size = size_requested

    size_x = size
    size_b = size + 1

    control = QuantumRegister(1, "c")
    reg_x = QuantumRegister(size_x, "x")
    reg_b = QuantumRegister(size_b, "b")
    # The one ancillary register for addition
    reg_anc = QuantumRegister(1, "anc")

    u = QuantumCircuit(control, reg_x, reg_b, reg_anc, name="CMUL(%d)MOD(%d)" % (a, N))

    # mult -> ax mod N
    u.append(mult_mod_N_partial_c(a, N, size_x, size_b), u.qubits)

    # controlled swaps
    # swap the bottom with the x register
    for i in range(0, size_x):
        u.cswap(control[0], reg_x[i], reg_b[i])

    # the last bit in b does not have to be swapped since it stays 0 after each addition
    inv_a = pow(a, -1, N)
    u.append(mult_mod_N_partial_c(inv_a, N, size_x, size_b).inverse(), u.qubits)

    return u


def mult_mod_N_partial_c(y, N, size_x, size_b):
    """Controlled out-of-place multiplication as described in
    Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.

    Register sizes:
    1 control
    size_x input register |x> -> |x>
    size_b result register |y> -> |y + a*x mod N>
    1 ancilla
    """
    control = QuantumRegister(1, "c")
    reg_x = QuantumRegister(size_x, "x")
    reg_b = QuantumRegister(size_b, "b")
    reg_anc = QuantumRegister(1, "anc")

    mult = QuantumCircuit(control, reg_x, reg_b, reg_anc, name="CMUL0(%d)MOD(%d)" % (y, N))

    # QFT on b
    mult.append(QFT(size_b, do_swaps=False), reg_b)

    # For each x-bit add a controlled addition with 2^k*y
    for k in range(0, size_x):
        adder_gate = phi_add_mod_N_cc(size_b, (2 ** k * y % N), N)
        # the gate at x_k is controlled by the global control and x_k
        # it works on register b + the ancillary bit
        controls = [control[0], reg_x[k]]
        mult.append(adder_gate, controls + list(reg_b) + [reg_anc[0]])

    # IQFT on b
    mult.append(QFT(size_b, do_swaps=False).inverse(), reg_b)

    return mult
