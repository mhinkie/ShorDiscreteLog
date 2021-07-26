from qiskit import QuantumCircuit, QuantumRegister
from .hrs_adder import c_adder_hrs2
from .hrs_carry import c_carry_gate, cc_carry_gate
import numpy as np


def c_comp(n, value):
    """Controlled comparator (with constant) using the controlled carry gate.

    Register sizes:
    1 control
    n input
    n-1 borrowed qubits
    1 result (= msb of subtraction)
    """
    control = QuantumRegister(1, "c")
    inreg = QuantumRegister(n, "in")
    breg = QuantumRegister(n - 1, "b")
    res = QuantumRegister(1, "res")

    comp = QuantumCircuit(control, inreg, breg, res, name="C_COMP(%d)" % value)

    valuesub = 0 - value

    c_bin = np.binary_repr(valuesub, n + 1)
    c_bin_rev = c_bin[::-1]  # reverse to have lsb first

    comp.append(c_carry_gate(n + 1, c_bin_rev), list(control) + list(inreg) + list(res) + list(breg))

    return comp


def cc_comp(n, value):
    """Double controlled comparator (with constant) using the double controlled carry gate.

    Register sizes:
    2 control
    n input
    n-1 borrowed qubits
    1 result (= msb of subtraction)
    """
    control = QuantumRegister(2, "c")
    inreg = QuantumRegister(n, "in")
    breg = QuantumRegister(n - 1, "b")
    res = QuantumRegister(1, "res")

    comp = QuantumCircuit(control, inreg, breg, res, name="CC_COMP(%d)" % value)

    valuesub = 0 - value
    # print("valuesub: ", valuesub)

    c_bin = np.binary_repr(valuesub, n + 1)
    c_bin_rev = c_bin[::-1]  # reverse to have lsb first
    # print("compare val: ", c_bin)

    comp.append(cc_carry_gate(n + 1, c_bin_rev), list(control) + list(inreg) + list(res) + list(breg))

    return comp


def c_carry_add_mod_hrs(n, a, N):
    """Modular addition (one control qubit) with a constant "a" as described in:
    Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication.

    Register sizes:
    1 control
    n input = output register.
    n-1 borrowed qubits for addition
    1 ancilla qubit for modular addition."""

    control = QuantumRegister(1, "c")
    breg = QuantumRegister(n, "b")
    greg = QuantumRegister(n - 1, "g")
    msbreg = QuantumRegister(1, "msb")  # msb for compare

    adder = QuantumCircuit(control, breg, greg, msbreg, name="CADD(%d)MOD(%d)" % (a, N))

    # CMP
    adder.append(c_comp(n, N - a), [control[0]] + list(breg) + list(greg) + [msbreg[0]])

    # Adder conditioned on msb according to cmp
    adder.append(c_adder_hrs2(n, a), [msbreg[0]] + list(breg) + [greg[0], greg[1]])  # can just borrow one of breg

    # flip msb according to control
    adder.cx(control, msbreg)

    # subtr
    adder.append(c_adder_hrs2(n, (a - N)), [msbreg[0]] + list(breg) + [greg[0], greg[1]])

    # If we added a before (without decreasing by N), then
    # it is >= a in any case: r >= a (and have msb = 1)
    # a compare flips the result bit if r < a
    # so if we flip the result bit before the last compare and then do the OP
    # it sets the result bit correct (this is why there is no conditional
    # flip to flip the bit back to what it was before
    adder.append(c_comp(n, a), [control[0]] + list(breg) + list(greg) + [msbreg[0]])

    return adder


def cc_carry_add_mod_hrs(n, a, N):
    """Modular addition (two control qubits) with a constant "a" as described in:
    Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication.

    Register sizes:
    2 control
    n input = output register.
    n-1 borrowed qubits for addition
    1 ancilla qubit for modular addition."""

    control = QuantumRegister(2, "c")
    breg = QuantumRegister(n, "b")
    greg = QuantumRegister(n - 1, "g")
    msbreg = QuantumRegister(1, "msb")  # msb for compare

    adder = QuantumCircuit(control, breg, greg, msbreg, name="CCADD(%d)MOD(%d)" % (a, N))

    # CMP
    adder.append(cc_comp(n, N - a), list(control) + list(breg) + list(greg) + [msbreg[0]])

    # Adder conditioned on msb according to cmp
    adder.append(c_adder_hrs2(n, a), [msbreg[0]] + list(breg) + [greg[0], greg[1]])  # can just borrow one of breg

    # flip msb according to control
    adder.ccx(control[0], control[1], msbreg)

    # subtr
    adder.append(c_adder_hrs2(n, (a - N)), [msbreg[0]] + list(breg) + [greg[0], greg[1]])

    # If we added a before (without decreasing by N), then
    # it is >= a in any case: r >= a (and have msb = 1)
    # a compare flips the result bit if r < a
    # so if we flip the result bit before the last compare and then do the OP
    # it sets the result bit correct (this is why there is no conditional
    # flip to flip the bit back to what it was before
    adder.append(cc_comp(n, a), list(control) + list(breg) + list(greg) + [msbreg[0]])

    return adder
