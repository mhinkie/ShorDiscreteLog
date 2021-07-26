from qiskit.circuit.library import QFT
from qiskit import QuantumCircuit, QuantumRegister
from .draper_add import phi_addition_cc, phi_addition_controlled_ignorebits, phi_addition_controlled


def mult_montgomery_partial_c(n, classical_y_montg, p):
    """Out of place montgomery multiplication using fouier space adder:
    |x>|0...0> -> |x>|0>|xy mod p>|0...0>

    Built after chapter 6.2 of:
    Rich Rines and Isaac Chuang. High performance quantum modular multipliers.

    The input for this stage (the constant y) is expected to be supplied
    in montgomery form (yR mod p).

    Internally performs: x * (yR mod p) and then performs montgomery reduction
    to xy mod p.

    Register sizes:
    1 control
    n input register |x>
    2n + 1 result + ancilla register (note: the result is not aligned to a border of the register
    for the out-of-place multiplication - this is reflected in the structure of the in-place version)"""
    ctrl = QuantumRegister(1, "ctrl")
    smallreg = QuantumRegister(n, "smallreg")
    bigreg = QuantumRegister(n + n + 1, "bigreg")  # + ancilla for reduction

    mult = QuantumCircuit(ctrl, smallreg, bigreg, name="CCPHIMULT_partial(%d mod %d)" % (classical_y_montg, p))
    mult.append(QFT(n + n + 1, do_swaps=False), bigreg)

    # Perform multiplication by repeated addition (controlled)
    for i in range(0, n):  # for each bit in x
        # controlled addition by 2^i Y' mod N
        summand = (2 ** i * classical_y_montg % p)
        mult.append(phi_addition_cc(n + n + 1, summand, 1), [ctrl[0], smallreg[i]] + list(bigreg))

    # now either |x>|t> or |x>|0>

    # montgomery reduction to get proper result
    # in each iteration the lsb will be fixed to the |u> register
    # it is also used as the control for the subtraction of N (or N-1 later)
    # the division by 2 is done implicitly as in each iteration a smaller register is used
    # lsb -> msb
    subtract_by = int(p)

    for i in range(0, n):
        # use the lsb as input
        qubits = [bigreg[i]] + list(bigreg[i + 1:])
        # print(qubits)
        # hadamard transform on the control bit to get it to computational basis
        mult.h(bigreg[i])

        # ignore the subtraction on bits that are already fixed in u
        sub_gate = phi_addition_controlled_ignorebits(n + n + 1 - i, subtract_by, -1, 1)  # always ignore the lsb
        mult.append(sub_gate, qubits)  # = subtraction because factor = -1

    mult.barrier()

    # transform back from qft space to extract sign
    mult.append(QFT(n + 1, do_swaps=False).inverse(), bigreg[n:])

    signbit = bigreg[-1:]

    # return to QFT for the remaining l bits
    mult.append(QFT(n, do_swaps=False), bigreg[n:n + n])

    # controlled addition of p (if sign is negative)
    mult.append(phi_addition_controlled(n, p, 1), [signbit] + list(bigreg[n:n + n]))

    # now the sign is flipped according to lsb of the result as before
    # as in QFTy the lsb is readable witout transformation, cnot suffices
    # z-variant: Have to apply hadamard on lsb, and then hadamard on lsb again
    mult.h(bigreg[n])
    mult.cx(bigreg[n], signbit)
    mult.h(bigreg[n])
    # the signbit will now be regarded as part of u again

    # now have to uncompute u' as before
    # the paper now gives a way to uncompute this u' using an equality (= subtractions)
    uncom_qubits = list(bigreg[:n]) + (bigreg[-1:])
    mult.append(QFT(n + 1, do_swaps=False), uncom_qubits)

    pinvbig = pow(p, -1, 2 ** (n + 1))  # p inverse mod 2^(l+1)
    for i in range(0, n):
        subtract_by = ((2 ** i * classical_y_montg % p) * pinvbig) % 2 ** (n + 1)
        mult.append(phi_addition_cc(n + 1, subtract_by, -1), [ctrl[0], smallreg[i]] + list(uncom_qubits))

    # TODO: i.m.o. hadamard would be enough for this register

    # qft for result register
    mult.append(QFT(n, do_swaps=False).inverse(), bigreg[n:n + n])

    # qft for ancillas
    mult.append(QFT(n + 1, do_swaps=False).inverse(), uncom_qubits)

    return mult


def mult_montgomery_c(n, y, p):
    """Controlled in-place montgomery multiplication
    |x>|0...0> -> |xy mod p>|0...0>

    Using two out-of-place multiplications.
    The input y is expected as usual (not in Montgomery form).

    Built after chapter 6.2 of:
    Rich Rines and Isaac Chuang. High performance quantum modular multipliers.

    Register sizes:
    1 control
    n input = output |x> -> |xy mod p>
    2n + 1 ancilla
    """
    ctrl = QuantumRegister(1, "ctrl")
    smallreg = QuantumRegister(n, "smallreg")
    bigreg = QuantumRegister(n + n + 1, "bigreg")  # + ancilla for reduction

    mult = QuantumCircuit(ctrl, smallreg, bigreg, name="CPHIMULT(%d mod %d)" % (y, p))

    R = 2 ** n
    # montgomery repr
    classical_y_montg = y * R % p

    # = controlled: |x>|0> -> |x>|0>|phi(xy mod p)>|0...0>
    mult.append(mult_montgomery_partial_c(n, classical_y_montg, p), mult.qubits)

    # swap result and input
    for i in range(0, n):
        mult.cswap(ctrl[0], smallreg[i], bigreg[n + i])
    # = |xy mod p>|0>|x>|0...0>

    inverse_y = pow(y, -1, p)  # should always work as p is prime
    inverse_y_montg = inverse_y * R % p

    # do the inverse of the mult with input the inverse of Y
    mult.append(mult_montgomery_partial_c(n, inverse_y_montg, p).inverse(), mult.qubits)
    # = |xy mod p>|0>|phi(0)>|0...0> = |xy mod p>|0>

    return mult
