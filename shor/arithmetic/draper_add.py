from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import QFT
import numpy as np


def phi_addition_cc(n, b, factor):
    """Double controlled addition in fourier space:
    |phi(x)> -> |phi(x + factor*b)>"""
    ctrl = QuantumRegister(2, "c")
    reg = QuantumRegister(n, "r")
    addc = QuantumCircuit(ctrl, reg, name="CCPHIADD(%d*%d)" % (b, factor))

    for k in range(0, n):
        addc.append(ccp(factor * get_angle(k, b)), list(ctrl) + [reg[k]])

    return addc



def ccp(theta):
    """Double controlled phase gates using single controlled phase gates."""

    ctrl = QuantumRegister(2, "ctrl")
    q = QuantumRegister(1, "q")
    ccp = QuantumCircuit(ctrl, q, name="CCP(%d)" % theta)

    ccp.cp(theta / 2, ctrl[1], q)
    ccp.cx(ctrl[0], ctrl[1])
    ccp.cp(-theta / 2, ctrl[1], q)
    ccp.cx(ctrl[0], ctrl[1])
    ccp.cp(theta / 2, ctrl[0], q)

    return ccp


def phi_addition_controlled(n, b, factor):
    """
    Single controlled addition in Fourier space:
    |phi(x)> -> |phi(x + factor*b)>

    See: Thomas G. Draper. Addition on a quantum computer.
    """
    return phi_addition_controlled_ignorebits(n, b, factor, 0)


def phi_addition_controlled_ignorebits(n, b, factor, ignorebits):
    """Single controlled addition in Fourier space:
    |phi(x)> -> |phi(x + factor*b)>

    ignorebits allows to remove all operations on the qubits < ignorebits
    (all less significant ones). Is used in the implementation of montgomery
    multiplication.

    Using approach from:
    Thomas G. Draper. Addition on a quantum computer.

    """
    ctrl = QuantumRegister(1, "c")
    reg = QuantumRegister(n - ignorebits, "r")
    addc = QuantumCircuit(ctrl, reg, name="CPHIADD(%d*%d,ignore=<%d)" % (b, factor, ignorebits))
    if ignorebits <= 0:
        addc = QuantumCircuit(ctrl, reg, name="CPHIADD(%d*%d)" % (b, factor))

        # addition in y-fourier space
    for k in range(0, n):
        if k < ignorebits:
            continue
        # rotate bit k by b*pi/(2^k+1)
        addc.cp(factor * get_angle(k, b), ctrl, reg[k - ignorebits])

    return addc



def add_draper(n, b, factor):
    """Adder |x> -> |x + factor * b> using adder in fourier space.

    QFT is performed without swaps (not necessary since the adder
    can just be appended exactly as needed for the result."""
    areg = QuantumRegister(n, "a")

    adder = QuantumCircuit(areg, name="ADD(%d*%d)" % (b, factor))

    # dont do swaps so QFT is as in paper
    adder.append(QFT(n, do_swaps=False), areg)

    adder.append(phi_add_builtin(n, b, factor), areg)

    adder.append(QFT(n, do_swaps=False).inverse(), areg)

    return adder


def get_angle(k, b):
    """Computes the angle of rotation based on the input constant.
    = b*pi/2^k"""
    # rotate by b*pi/2^k
    return b * np.pi / (2 ** k)


def phi_add_builtin(n, b, factor):
    """Adder |phi(x)> -> |phi(x + factor*b)>.
    As opposed to phi_add_builtin_singlegates, only one rotation is performed
    for each qubit (since b is fully known in advance this rotation can be
    precomputed).

    See: Thomas G. Draper. Addition on a quantum computer."""
    areg = QuantumRegister(n, "phi_a")

    add = QuantumCircuit(areg, name="phiADD(%d*%d)" % (b, factor))

    for k in range(0, n):
        angle = get_angle(k, b)
        add.p(factor * angle, areg[k])

    return add


def phi_add_builtin_singlegates(n, b, factor):
    """
    |phi(x)> -> |phi(x + factor*b)>

    b and factor are builtin to the circuit.

    For debug/visualisation: Adder in Fourier space, where each
    input bit of the constant is treated as one qubit (although 'it is
    not there') - mimics the layout in:
    Thomas G. Draper. Addition on a quantum computer.

    factor allows simple creation of a subtraction (factor = -1)."""
    areg = QuantumRegister(n, "phi_a")

    add = QuantumCircuit(areg, name="phiADD(%d*%d)" % (b, factor))

    b_bin = np.binary_repr(b, n)
    # print("b: %s" % b_bin)
    # reverse the binary repr so the accessing the required values is easier
    b_bin = b_bin[::-1]

    for k in range(0, n):
        # apply all angles for a_k
        for i in range(0, k + 1):
            if b_bin[i] == '1':
                r = k - i + 1
                add.p(factor * 2 * np.pi / 2 ** r, areg[k])

    return add
