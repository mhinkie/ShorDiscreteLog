from qiskit import QuantumCircuit, QuantumRegister
from .incrementer import cinc
from .hrs_carry import c_carry_gate
import numpy as np
import math


def c_divide_and_conquer_adder(adder, n, xreg, o_borrowed, c_bin_rev, control_reg, borrowed_extra):
    """Controlled adder using divide and conquer approach and carry-gate as described in:
    Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication."""

    # with control bit
    # if there are only 2 qubits left (n=2), the carry from the high bit is
    # 1 if c_0 = 1 and the 0 otherwise => communicating the carry is just
    # flipping the low bit according to c_0 and the high bit
    # the borrowed bit does not have to be touched
    if n== 2:
        # print("In base case n=2")
        if c_bin_rev[0] == '1':
            # there is a carry if the low bit is 1
            adder.ccx(control_reg, xreg[0], xreg[1])
        else:
            pass  # there cannot be a carry if c_0 = 0

        # now just do the simple addition (= bit flip)
        # without carry (the carries are all accounted for)
        if c_bin_rev[0] == '1':
            adder.cx(control_reg, xreg[0])

        if c_bin_rev[1] == '1':
            adder.cx(control_reg, xreg[1])

        return

    # if this gets called with only one qubit (= if the previous had 3 for example)
    # the carry is already accounted for and we only need to flip according to the input
    if n == 1:
        # print("In base case n=1")
        if c_bin_rev[0] == '1':
            adder.cx(control_reg, xreg[0])

        return

    nlow = math.ceil((n) / 2)
    nhigh = n - nlow

    # print("nlow: ", nlow)
    # print("nhigh: ", nhigh)

    # to be able to handle the toggle
    # increment high register before
    # incrementer for nhigh qubits with 1 control bit needs nhigh+1 qubits borrowed
    # in the case where n/2 is integer, nhigh = nlow, and there are not enough borrowed qubits
    # to perform increment - here one extra borrowed qubit is used (from outside)
    num_borrowed_needed = nhigh + 1
    borrowed_for_incrementer = xreg[:nlow]
    if num_borrowed_needed > len(borrowed_for_incrementer):
        borrowed_for_incrementer.append(borrowed_extra)
    # print([o_borrowed]+xreg[nlow:]+borrowed_for_incrementer)
    adder.append(cinc(nhigh), [o_borrowed] + xreg[nlow:] + borrowed_for_incrementer)

    # do it for each bit seperately
    for i, bit in enumerate(xreg[nlow:]):
        adder.cx(o_borrowed, bit)

    # Add carry computation
    # compute_highest_bit_sum has one register for the input, where the highest of those is for the carry
    # print("appending low bit adder for: %s (reversed)" % c_bin_rev[:nlow])
    # 0 char is added in back because compute_highest_bit_sum expects some bit for carry bit aswell
    # for example - perform the addition for a 4-bit number although I only have 3 bits - the 4th bit
    # is then purely the carry
    num_borrowed = nlow - 1  # number of borrowed qubits for carry
    # print("need to borrow: ", num_borrowed)
    adder.append(c_carry_gate(nlow + 1, c_bin_rev[:nlow] + '0'),
                 [control_reg] + list(xreg[:nlow]) + [o_borrowed] + list(xreg[nlow:nlow + num_borrowed]))

    # o_borrowed is now toggled if there is a carry

    adder.append(cinc(nhigh), [o_borrowed] + xreg[nlow:] + borrowed_for_incrementer)

    adder.append(c_carry_gate(nlow + 1, c_bin_rev[:nlow] + '0'),
                 [control_reg] + list(xreg[:nlow]) + [o_borrowed] + list(xreg[nlow:nlow + num_borrowed]))

    for i, bit in enumerate(xreg[nlow:]):
        adder.cx(o_borrowed, bit)

    # now recursively do the same for x_L and x_H
    # print("calling recursive with size %d for %s" % (nlow, xreg[:nlow]))
    c_divide_and_conquer_adder(adder, nlow, xreg[:nlow], o_borrowed, c_bin_rev[:nlow], control_reg, borrowed_extra)
    # print("calling recursive with size %d for %s" % (nhigh, xreg[nlow:nlow+nhigh]))
    c_divide_and_conquer_adder(adder, nhigh, xreg[nlow:nlow + nhigh], o_borrowed, c_bin_rev[nlow:nlow + nhigh],
                               control_reg, borrowed_extra)


def c_adder_hrs2(n, c):
    """Controlled adder using the divide and conquer adder."""
    control = QuantumRegister(1, "c")
    xreg = QuantumRegister(n, "x")
    o_borrowed = QuantumRegister(2, "o_b")  # one extra borrowed qubit for incrementer = 2

    # if a negative number is converted using binary_repr, numpy might
    # complain that the supplied size (n) is too small, but this is
    # expected. To get rid of the warning the number is converted with one
    # additional bit, and the first one (the extra msb) is then ignored
    c_bin = np.binary_repr(c, n+1)
    c_bin = c_bin[1:]
    # print("C in binary: ", c_bin)
    c_bin_rev = c_bin[::-1]  # reverse to have lsb first

    adder = QuantumCircuit(control, xreg, o_borrowed, name="ADD(%d)" % c)

    c_divide_and_conquer_adder(adder, n, xreg, o_borrowed[0], c_bin_rev, control, o_borrowed[1])

    return adder

