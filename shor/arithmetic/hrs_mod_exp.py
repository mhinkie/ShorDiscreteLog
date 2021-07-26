from qiskit import QuantumCircuit, QuantumRegister
from .hrs_mod_mult import c_mult_hrs


def mod_exp_hrs(ntop, nbot, a, N):
    """Modular exponentiation using controlled multipliers.
    |x>|y>|0...0> -> |x>|y * a^x mod N>|0...0>

    Register sizes:
    ntop top register |x> -> |x>
    nbot bottom register |y> -> |y*a^x mod N>
    nbot+1 ancilla |0> -> |0>

    See: Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication.
    and
    Stephane Beauregard. Circuit for Shor’s algorithm using 2n+3 qubits.
    """

    top = QuantumRegister(ntop, "top")
    bot = QuantumRegister(2 * nbot + 1, "bot")

    circ = QuantumCircuit(top, bot, name="%d^x MOD(%d)" % (a, N))

    for i in range(0, ntop):
        # top[i] is the control bit
        circ.append(c_mult_hrs(nbot, a ** (2 ** i) % N, N), [top[i]] + list(bot))

    return circ
