from qiskit import QuantumCircuit, QuantumRegister


def cinc(n):
    """Incrementer using borrowed qubits as described in:
    Thomas Häner, Martin Roetteler, and Krysta M. Svore. Factoring using 2n+2 qubits
    with Toffoli based modular multiplication."""
    ctrl = QuantumRegister(1, "control")
    xreg = QuantumRegister(n, "x")
    breg = QuantumRegister(n + 1, "borrowed")

    inc = QuantumCircuit(ctrl, xreg, breg, name="C+1")

    subtractor = adder_ttk(n + 1,
                           handle_overflow=False).inverse()  # |a>|a+c> -> |a>|c> (first register n qubits, second n+1)

    logical_xreg = [ctrl] + list(xreg)

    # |x>|g> -> |x-g>|g>
    inc.append(subtractor, list(breg) + list(logical_xreg))

    # not on g register gives g'-1, where g' is the twos complement g
    inc.x(breg)

    # subtract again |x-g-g'+1>|g'-1>
    inc.append(subtractor, list(breg) + list(logical_xreg))

    # invert borrowed qubits again
    inc.x(breg)

    # x on control
    inc.x(ctrl)

    return inc


def adder_ttk(n, alternating_reg=False, handle_overflow=True):
    """|a,b> -> |a, a+b>, where b has one extra qubit.
    Based on: Yasuhiro Takahashi, Seiichiro Tani, and Noboru Kunihiro.
    Quantum addition circuits and unbounded fan-out.
    """
    areg = QuantumRegister(n, "a")
    if handle_overflow:  # ignores the last qubit (= carry after addition)
        breg = QuantumRegister(n + 1, "b")
    else:
        breg = QuantumRegister(n, "b")

    # only for visualization: alternate registers so same i is together in output
    if (alternating_reg):
        qubits = []
        for i in range(0, n):
            qubits.append(breg[i])
            qubits.append(areg[i])

        if handle_overflow:
            qubits.append(breg[n])
        add = QuantumCircuit(qubits, name="ADD(q,q)")
    else:
        add = QuantumCircuit(areg, breg, name="ADD(q,q)")

    for i in range(1, n):
        add.cx(areg[i], breg[i])

    if handle_overflow:
        add.cx(areg[n - 1], breg[n])

    for i in reversed(range(1, n - 1)):
        add.cx(areg[i], areg[i + 1])

    for i in range(0, n - 1):
        add.ccx(breg[i], areg[i], areg[i + 1])

    if handle_overflow:
        add.ccx(breg[n - 1], areg[n - 1], breg[n])

    for i in reversed(range(1, n)):
        add.cx(areg[i], breg[i])
        add.ccx(breg[i - 1], areg[i - 1], areg[i])

    for i in range(1, n - 1):
        add.cx(areg[i], areg[i + 1])

    for i in range(0, n):
        add.cx(areg[i], breg[i])

    return add
