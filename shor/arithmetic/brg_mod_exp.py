from qiskit import QuantumCircuit, QuantumRegister
import math
from .brg_mod_mult import mult_mod_N_c


def mod_exp_brg(n, a, N):
    """Modular exponentiation:
    |x>|y>|0...0> -> |x>|y * a^x mod N>|0...0>

    Register sizes:
    n top register |x> -> |x>
    n bottom register |y> -> |y*a^x mod N>
    n+2 ancilla |0> -> |0>

    See: Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.
    """

    # mult register has size 2n+3, where 1 is the control bit
    size_input = 2*math.ceil(math.log(N, 2))+3

    x_reg = QuantumRegister(n, "x")
    input_reg = QuantumRegister(size_input-1, "input") # one bit of the u-gate is control

    circ = QuantumCircuit(x_reg, input_reg, name="%d^x MOD(%d)" % (a, N))

    for i in range(0, n):
        mulgate = mult_mod_N_c(a**(2**i) % N, N)
        circ.append(mulgate, [x_reg[i]] + list(input_reg))

    return circ
