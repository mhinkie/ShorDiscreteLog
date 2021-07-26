from qiskit import QuantumCircuit, QuantumRegister

"""carry gate calculates the carry after adding a n-bit constant to a n-bit
register.
The n-bit constant is expected as a reversed binary string
(reversed because then the lsb is first and can be indexed with [0])"""


def carry_gate(n, c_bin_rev):
    areg = QuantumRegister(n, "a")
    greg_actual = QuantumRegister(n - 2, "g")

    carryc = QuantumCircuit(areg, greg_actual, name="calc_msb(n=%d,c=%s)" % (n, c_bin_rev[::-1]))

    greg = [areg[0]] + list(greg_actual)  # optimization: the first carry is actually just the areg

    c_bin = c_bin_rev

    #################### FIRST FORWARD ##############

    # I start with the last bit (= result bit)
    # since this bit is assumed to be 0 (since I use one more bit than necessary),
    carryc.cx(greg[n - 2], areg[n - 1])

    # Since we apply the controlled operations, before the flips are actually done, we start with the highest bit
    # (see above) - the last one and first one will get extra treatment
    for i in reversed(range(1, n - 1)):
        if i == 1:
            # if a0 = 0 the carry from 1 to 2 is influenced only
            # by a1 (since there is no carry from 0 to 1)
            # otherwise proceed as usual
            if c_bin[0] == '0':
                if c_bin[i] == '1':
                    carryc.cx(areg[i], greg[i])
                    carryc.x(areg[i])
                    # the ccx gate is removed (would cancel out with down dir)
                else:
                    pass  # in the case that this is 0 aswell, there also will
                    # be no ccx
                continue

        if c_bin[i] == '1':
            # if c is 1, we have cnot and x
            carryc.cx(areg[i], greg[i])
            carryc.x(areg[i])
            carryc.ccx(greg[i - 1], areg[i], greg[i])
        else:
            # if c is 0, we only have one ccnot gate
            carryc.ccx(greg[i - 1], areg[i], greg[i])

    # for the first bit we can assume that the incoming carry is 0
    # since the first carry is only influenced by the first bit,
    # since greg actually includes the first a-register this is already respected
    # so for a0 nothing is done

    # now move from lsb to msb and repeat the operations that are controlled by the carries.
    # Now we can be sure that they are only applied if they were actually toggled
    for i in range(1, n - 1):
        # since greg[0] is actually not a carry anymore, we dont have to apply it for this
        if (i - 1) == 0:
            continue
        # in both cases c_i = 1 and c_i = 0, only the ccnot gate is controlled by g
        # therefore I only need to apply this again
        carryc.ccx(greg[i - 1], areg[i], greg[i])

    # and as before the last bit is applied according to the toggle g
    carryc.cx(greg[n - 2], areg[n - 1])

    # in preperation for the rest of the algorithm, I already add the x-gate to the last bit
    # if msb of c is 1 (this is the same as just adding the carry to msb)
    if c_bin[n - 1] == '1':
        carryc.x(areg[n - 1])

    #################### BACKWARD, without touching msb of a #############

    # reverse ccnot from msb to lsb
    for i in reversed(range(1, n - 1)):
        # as before greg[0] is no actual carry register
        # => dont apply if greg[i-1] = greg[0]
        if (i - 1) == 0:
            continue
        carryc.ccx(greg[i - 1], areg[i], greg[i])

    # not needed due to optimization
    # reverse op on first carry
    # if c_bin[0] == '1':
    #    carryc.cx(areg[0], greg[0])

    # reverse ccnot + cnot + x from lsb to msb
    for i in range(1, n - 1):
        if i == 1:
            # if a0 = 0 the carry from 1 to 2 is influenced only
            # by a1 (since there is no carry from 0 to 1)
            # otherwise proceed as usual
            if c_bin[0] == '0':
                if c_bin[i] == '1':
                    carryc.x(areg[i])
                    carryc.cx(areg[i], greg[i])
                    # the ccx gate is removed (would cancel out with down dir)
                else:
                    pass  # in the case that this is 0 aswell, there also will
                    # be no ccx
                continue

        if c_bin[i] == '1':
            # if c is 1, we have cnot and x (reverse order)
            carryc.ccx(greg[i - 1], areg[i], greg[i])
            carryc.x(areg[i])
            carryc.cx(areg[i], greg[i])
        else:
            # if c is 0, we only have one ccnot gate
            carryc.ccx(greg[i - 1], areg[i], greg[i])

    return carryc


def c_carry_gate(n, c_bin_rev):
    # Controlled carry gate (controlled by first bit)
    # control all operations on the output bit a[n-1]
    c = QuantumRegister(1, "control")
    areg = QuantumRegister(n, "a")
    greg_actual = QuantumRegister(n - 2, "g")

    carryc = QuantumCircuit(c, areg, greg_actual, name="calc_msb(n=%d,c=%s)" % (n, c_bin_rev[::-1]))

    greg = [areg[0]] + list(greg_actual)  # optimization: the first carry is actually just the areg

    c_bin = c_bin_rev

    #################### FIRST FORWARD ##############

    # I start with the last bit (= result bit)
    # since this bit is assumed to be 0 (since I use one more bit than necessary),
    carryc.ccx(c, greg[n - 2], areg[n - 1])

    # Since we apply the controlled operations, before the flips are actually done, we start with the highest bit
    # (see above) - the last one and first one will get extra treatment
    for i in reversed(range(1, n - 1)):
        if i == 1:
            # if a0 = 0 the carry from 1 to 2 is influenced only
            # by a1 (since there is no carry from 0 to 1)
            # otherwise proceed as usual
            if c_bin[0] == '0':
                if c_bin[i] == '1':
                    carryc.cx(areg[i], greg[i])
                    carryc.x(areg[i])
                    # the ccx gate is removed (would cancel out with down dir)
                else:
                    pass  # in the case that this is 0 aswell, there also will
                    # be no ccx
                continue

        if c_bin[i] == '1':
            # if c is 1, we have cnot and x
            carryc.cx(areg[i], greg[i])
            carryc.x(areg[i])
            carryc.ccx(greg[i - 1], areg[i], greg[i])
        else:
            # if c is 0, we only have one ccnot gate
            carryc.ccx(greg[i - 1], areg[i], greg[i])

    # for the first bit we can assume that the incoming carry is 0
    # since the first carry is only influenced by the first bit,
    # since greg actually includes the first a-register this is already respected
    # so for a0 nothing is done

    # now move from lsb to msb and repeat the operations that are controlled by the carries.
    # Now we can be sure that they are only applied if they were actually toggled
    for i in range(1, n - 1):
        # since greg[0] is actually not a carry anymore, we dont have to apply it for this
        if (i - 1) == 0:
            continue
        # in both cases c_i = 1 and c_i = 0, only the ccnot gate is controlled by g
        # therefore I only need to apply this again
        carryc.ccx(greg[i - 1], areg[i], greg[i])

    # and as before the last bit is applied according to the toggle g
    carryc.ccx(c, greg[n - 2], areg[n - 1])

    # in preperation for the rest of the algorithm, I already add the x-gate to the last bit
    # if msb of c is 1 (this is the same as just adding the carry to msb)
    if c_bin[n - 1] == '1':
        carryc.cx(c, areg[n - 1])

    #################### BACKWARD, without touching msb of a #############

    # reverse ccnot from msb to lsb
    for i in reversed(range(1, n - 1)):
        # as before greg[0] is no actual carry register
        # => dont apply if greg[i-1] = greg[0]
        if (i - 1) == 0:
            continue
        carryc.ccx(greg[i - 1], areg[i], greg[i])

    # not needed due to optimization
    # reverse op on first carry
    # if c_bin[0] == '1':
    #    carryc.cx(areg[0], greg[0])

    # reverse ccnot + cnot + x from lsb to msb
    for i in range(1, n - 1):
        if i == 1:
            # if a0 = 0 the carry from 1 to 2 is influenced only
            # by a1 (since there is no carry from 0 to 1)
            # otherwise proceed as usual
            if c_bin[0] == '0':
                if c_bin[i] == '1':
                    carryc.x(areg[i])
                    carryc.cx(areg[i], greg[i])
                    # the ccx gate is removed (would cancel out with down dir)
                else:
                    pass  # in the case that this is 0 aswell, there also will
                    # be no ccx
                continue

        if c_bin[i] == '1':
            # if c is 1, we have cnot and x (reverse order)
            carryc.ccx(greg[i - 1], areg[i], greg[i])
            carryc.x(areg[i])
            carryc.cx(areg[i], greg[i])
        else:
            # if c is 0, we only have one ccnot gate
            carryc.ccx(greg[i - 1], areg[i], greg[i])

    return carryc


def cc_carry_gate(n, c_bin_rev):
    # Controlled carry gate (controlled by first bit)
    # control all operations on the output bit a[n-1]
    c = QuantumRegister(2, "control")
    areg = QuantumRegister(n, "a")
    greg_actual = QuantumRegister(n - 2, "g")

    carryc = QuantumCircuit(c, areg, greg_actual, name="calc_msb(n=%d,c=%s)" % (n, c_bin_rev[::-1]))

    greg = [areg[0]] + list(greg_actual)  # optimization: the first carry is actually just the areg

    c_bin = c_bin_rev

    #################### FIRST FORWARD ##############

    # I start with the last bit (= result bit)
    # since this bit is assumed to be 0 (since I use one more bit than necessary),
    carryc.mct(list(c) + [greg[n - 2]], areg[n - 1])

    # Since we apply the controlled operations, before the flips are actually done, we start with the highest bit
    # (see above) - the last one and first one will get extra treatment
    for i in reversed(range(1, n - 1)):
        if i == 1:
            # if a0 = 0 the carry from 1 to 2 is influenced only
            # by a1 (since there is no carry from 0 to 1)
            # otherwise proceed as usual
            if c_bin[0] == '0':
                if c_bin[i] == '1':
                    carryc.cx(areg[i], greg[i])
                    carryc.x(areg[i])
                    # the ccx gate is removed (would cancel out with down dir)
                else:
                    pass  # in the case that this is 0 aswell, there also will
                    # be no ccx
                continue

        if c_bin[i] == '1':
            # if c is 1, we have cnot and x
            carryc.cx(areg[i], greg[i])
            carryc.x(areg[i])
            carryc.ccx(greg[i - 1], areg[i], greg[i])
        else:
            # if c is 0, we only have one ccnot gate
            carryc.ccx(greg[i - 1], areg[i], greg[i])

    # for the first bit we can assume that the incoming carry is 0
    # since the first carry is only influenced by the first bit,
    # since greg actually includes the first a-register this is already respected
    # so for a0 nothing is done

    # now move from lsb to msb and repeat the operations that are controlled by the carries.
    # Now we can be sure that they are only applied if they were actually toggled
    for i in range(1, n - 1):
        # since greg[0] is actually not a carry anymore, we dont have to apply it for this
        if (i - 1) == 0:
            continue
        # in both cases c_i = 1 and c_i = 0, only the ccnot gate is controlled by g
        # therefore I only need to apply this again
        carryc.ccx(greg[i - 1], areg[i], greg[i])

    # and as before the last bit is applied according to the toggle g
    carryc.mct(list(c) + [greg[n - 2]], areg[n - 1])

    # in preperation for the rest of the algorithm, I already add the x-gate to the last bit
    # if msb of c is 1 (this is the same as just adding the carry to msb)
    if c_bin[n - 1] == '1':
        carryc.ccx(c[0], c[1], areg[n - 1])

    #################### BACKWARD, without touching msb of a #############

    # reverse ccnot from msb to lsb
    for i in reversed(range(1, n - 1)):
        # as before greg[0] is no actual carry register
        # => dont apply if greg[i-1] = greg[0]
        if (i - 1) == 0:
            continue
        carryc.ccx(greg[i - 1], areg[i], greg[i])

    # not needed due to optimization
    # reverse op on first carry
    # if c_bin[0] == '1':
    #    carryc.cx(areg[0], greg[0])

    # reverse ccnot + cnot + x from lsb to msb
    for i in range(1, n - 1):
        if i == 1:
            # if a0 = 0 the carry from 1 to 2 is influenced only
            # by a1 (since there is no carry from 0 to 1)
            # otherwise proceed as usual
            if c_bin[0] == '0':
                if c_bin[i] == '1':
                    carryc.x(areg[i])
                    carryc.cx(areg[i], greg[i])
                    # the ccx gate is removed (would cancel out with down dir)
                else:
                    pass  # in the case that this is 0 aswell, there also will
                    # be no ccx
                continue

        if c_bin[i] == '1':
            # if c is 1, we have cnot and x (reverse order)
            carryc.ccx(greg[i - 1], areg[i], greg[i])
            carryc.x(areg[i])
            carryc.cx(areg[i], greg[i])
        else:
            # if c is 0, we only have one ccnot gate
            carryc.ccx(greg[i - 1], areg[i], greg[i])

    return carryc
