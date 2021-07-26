from .montgomery_mult import mult_montgomery_c
from qiskit import QuantumCircuit, QuantumRegister


def mod_exp_montgomery(n, a, p):
    """Modular exponentiation using montgomery multiplication.
    |x>|y>|0...0> -> |x>|y * a^x mod N>|0...0>

    See: Rich Rines and Isaac Chuang. High performance quantum modular multipliers.

    Register sizes:
    n top register |x> -> |x>
    n bottom register |y> -> |y*a^x mod N>
    2n+1 ancilla |0> -> |0>
    """

    top = QuantumRegister(n, "top")
    bot = QuantumRegister(3 * n + 1, "bot")

    circ = QuantumCircuit(top, bot, name="%d^x MOD(%d)" % (a, p))

    for i in range(0, n):
        # top[i] is the control bit
        circ.append(mult_montgomery_c(n, a ** (2 ** i) % p, p), [top[i]] + list(bot))

    return circ
