from typing import Dict, Optional, Union, Callable

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.aqua import QuantumAlgorithm, QuantumInstance
from qiskit.circuit.library import QFT
from qiskit.providers import Backend, BaseBackend

from shor.arithmetic.brg_mod_exp import mod_exp_brg

import math
import numpy as np
from abc import ABC, abstractmethod
from fractions import Fraction

from shor.arithmetic.brg_mod_mult import mult_mod_N_c


def mod_mul_brg(n, a, p):
    return mult_mod_N_c(a, p)


class DiscreteLogMoscaEkert(QuantumAlgorithm):
    def __init__(self,
                 b: int,
                 g: int,
                 p: int,
                 r: int = -1,
                 n: int = -1,
                 mod_exp_constructor: Callable[[int, int, int], QuantumCircuit] = None,
                 full_run: bool = False,
                 quantum_instance: Optional[
                     Union[QuantumInstance, Backend, BaseBackend, Backend]] = None) -> None:
        """
        Args:
            b: Finds discrete logarithm of b with respect to generator g and module p.
            g: Generator.
            p: Prime module.
            r: The order of g if it is known (otherwise it will be calculated)
            n: The size of the top register, if not given it will be inferred from the module p
            mod_exp_constructor: Function returns modular exponentiation circuit: (n, a, p) -> QuantumCircuit
                (n = size of top register)
            run_full (default: False), if set to True, the algorithm will not stop after finding the result,
                but will examine all results and log the probability
            quantum_instance: Quantum Instance or Backend
        """
        super().__init__(quantum_instance)
        self._b = b
        self._g = g
        self._p = p
        self._r = r
        self._n = n
        self._full_run = full_run

        if mod_exp_constructor is None:
            self._mod_exp_constructor = mod_exp_brg
        else:
            self._mod_exp_constructor = mod_exp_constructor

    @abstractmethod
    def _exec_qc(self, n, b, g, p, r, mod_exp_constructor, quantum_instance, shots) -> [(int, int, int)]:
        pass

    def get_circuit(self):
        if self._n == -1:
            n = math.ceil(math.log(self._p, 2))
        else:
            n = self._n
        print("constructing with size n=", n)
        return self._construct_circuit(n, self._b, self._g, self._p, self._r, self._mod_exp_constructor)

    @abstractmethod
    def _construct_circuit(self, n: int, b: int, g: int, p: int, r: int, mod_exp_constructor) -> QuantumCircuit:
        pass

    @staticmethod
    def find_period(result_list: [(int, int, int)], n: int, p: int, g: int) -> int:
        """The order of the generator is not given, it has to be determined.
        The list of results for the whole circuit will be used, since the
        measurement of stage1 allows for phase estimation of the operator
        (similar to Shor's algorithm for prime factorization)."""
        smallest_fitting_denominator = p+1  # init to p+1 (sth larger than the real result)
        for (y1, _, _) in result_list:
            meas_div = y1/(2**n)
            frac_meas = Fraction(meas_div).limit_denominator(p - 1)

            # check if denominator fits
            r_candidate = frac_meas.denominator
            if g**r_candidate % p == 1:
                # fits
                if r_candidate < smallest_fitting_denominator:
                    smallest_fitting_denominator = r_candidate

        print("Found r=", smallest_fitting_denominator)
        return smallest_fitting_denominator

    def _run(self) -> Dict:
        # construct circuit

        # size of top register is enough to hold all values mod p
        if self._n == -1:
            n = math.ceil(math.log(self._p, 2))
        else:
            n = self._n
        print("constructing with size n=", n)
        # print("top register size: ", n)

        if self._quantum_instance is not None:
            shots = 100

            result_list = self._exec_qc(n, self._b, self._g, self._p, self._r, self._mod_exp_constructor,
                                        self._quantum_instance, shots)

            num_fail = 0
            num_success = 0

            correct_m = -1

            for (y1, y2, freq) in result_list:
                if self._r == -1:
                    # period has to be calculated
                    self._r = self.find_period(result_list, n, self._p, self._g)

                # k is inferred from the measurement result from first stage
                k = int(round((y1 * self._r) / (2 ** n), 0))
                # print("k=", k)

                # m (the discrete log) is calculated by using the congruence:
                # ((km mod r)/r)*2^n = m_stage2
                # v = m_stage2*r/2^n
                v = (y2 * self._r) / (2 ** n)
                # print("v=", v)  # = km mod r

                # k inverse exists?
                if np.gcd(k, self._r) != 1:
                    # print("No inverse for k = fail")
                    num_fail += freq
                else:
                    kinv = pow(k, -1, self._r)
                    # print("kInv=", kinv)
                    m = int(round(v * kinv % self._r, 0))
                    # print("found m=", m)

                    # check if this m actually fits
                    if (self._g ** m % self._p) == self._b:
                        num_success += freq
                        # print("found correct m!")
                        correct_m = m
                        if not self._full_run:
                            break
                    else:
                        # print("found wrong m!")
                        num_fail += freq

            print("num_fail:", num_fail)
            print("num_success:", num_success)

            return {"m": correct_m, "success_prob": num_success / shots}

        return {"m": -1, "success_prob": 0}


class DiscreteLogMoscaEkertSharedRegister(DiscreteLogMoscaEkert):
    def __init__(self,
                 b: int,
                 g: int,
                 p: int,
                 r: int = -1,
                 n: int = -1,
                 mod_exp_constructor: Callable[[int, int, int], QuantumCircuit] = None,
                 full_run: bool = False,
                 quantum_instance: Optional[
                     Union[QuantumInstance, Backend, BaseBackend, Backend]] = None) -> None:
        super().__init__(b, g, p, r, n, mod_exp_constructor, full_run, quantum_instance)

    def _exec_qc(self, n, b, g, p, r, mod_exp_constructor, quantum_instance, shots) -> [(int, int, int)]:
        me_circuit = transpile(self._construct_circuit(n, b, g, p, r, mod_exp_constructor),
                               quantum_instance.backend)
        counts = quantum_instance.backend.run(me_circuit, shots=shots).result().get_counts(me_circuit)

        res = list()
        for result in counts.keys():
            # split result measurements
            result_s = result.split(" ")
            m_stage1 = int(result_s[1], 2)
            m_stage2 = int(result_s[0], 2)

            res.append((m_stage1, m_stage2, counts[result]))

        return res

    def _construct_circuit(self, n: int, b: int, g: int, p: int, r: int, mod_exp_constructor) -> QuantumCircuit:
        # infer size of circuit from modular exponentiation circuit
        mod_exp_g = mod_exp_constructor(n, g, p)
        mod_exp_b = mod_exp_constructor(n, b, p)

        iqft = QFT(n).inverse()

        total_circuit_qubits = mod_exp_g.num_qubits
        bottom_register_qubits = total_circuit_qubits - n

        topreg = QuantumRegister(n, "top")
        botreg = QuantumRegister(bottom_register_qubits, "bot")
        meas_stage1 = ClassicalRegister(n, "m1")
        meas_stage2 = ClassicalRegister(n, "m2")

        me_circuit = QuantumCircuit(topreg, botreg, meas_stage1, meas_stage2)

        # H on top
        me_circuit.h(topreg)

        # 1 on bottom
        me_circuit.x(botreg[0])

        # mod exp g^x mod p
        me_circuit.append(mod_exp_g, me_circuit.qubits)

        # iqft top
        me_circuit.append(iqft, topreg)

        # measure top register (stage 1)
        me_circuit.measure(topreg, meas_stage1)

        # reset top register
        me_circuit.reset(topreg)

        # h on top again
        me_circuit.h(topreg)

        # mod exp b^x' mod p
        me_circuit.append(mod_exp_b, me_circuit.qubits)

        # iqft top
        me_circuit.append(iqft, topreg)

        # measurement stage 2
        me_circuit.measure(topreg, meas_stage2)

        #print(me_circuit)

        return me_circuit


class DiscreteLogMoscaEkertSeperateRegister(DiscreteLogMoscaEkert):
    def __init__(self,
                 b: int,
                 g: int,
                 p: int,
                 r: int = -1,
                 n: int = -1,
                 mod_exp_constructor: Callable[[int, int, int], QuantumCircuit] = None,
                 full_run: bool = False,
                 quantum_instance: Optional[
                     Union[QuantumInstance, Backend, BaseBackend, Backend]] = None) -> None:
        super().__init__(b, g, p, r, n, mod_exp_constructor, full_run, quantum_instance)

    def _exec_qc(self, n, b, g, p, r, mod_exp_constructor, quantum_instance, shots) -> [(int, int, int)]:
        me_circuit = transpile(self._construct_circuit(n, b, g, p, r, mod_exp_constructor),
                               quantum_instance.backend)
        counts = quantum_instance.backend.run(me_circuit, shots=shots).result().get_counts(me_circuit)

        res = list()
        for result in counts.keys():
            # split result measurements
            result_s = result.split(" ")
            m_stage1 = int(result_s[1], 2)
            m_stage2 = int(result_s[0], 2)

            res.append((m_stage1, m_stage2, counts[result]))

        return res

    def _construct_circuit(self, n: int, b: int, g: int, p: int, r: int, mod_exp_constructor) -> QuantumCircuit:
        # infer size of circuit from modular exponentiation circuit
        mod_exp_g = mod_exp_constructor(n, g, p)
        mod_exp_b = mod_exp_constructor(n, b, p)

        iqft = QFT(n).inverse()

        total_circuit_qubits = mod_exp_g.num_qubits
        bottom_register_qubits = total_circuit_qubits - n

        top1reg = QuantumRegister(n, "topstage1")
        top2reg = QuantumRegister(n, "topstage2")
        botreg = QuantumRegister(bottom_register_qubits, "bot")
        meas_stage1 = ClassicalRegister(n, "m1")
        meas_stage2 = ClassicalRegister(n, "m2")

        me_circuit = QuantumCircuit(top1reg, top2reg, botreg, meas_stage1, meas_stage2)

        # H on top
        me_circuit.h(top1reg)

        # 1 on bottom
        me_circuit.x(botreg[0])

        # mod exp g^x mod p
        me_circuit.append(mod_exp_g, list(top1reg) + list(botreg))

        # iqft top
        me_circuit.append(iqft, top1reg)

        # h on top2
        me_circuit.h(top2reg)

        # mod exp b^x' mod p
        me_circuit.append(mod_exp_b, list(top2reg) + list(botreg))

        # iqft top
        me_circuit.append(iqft, top2reg)

        # measure top register (stage 1)
        me_circuit.measure(top1reg, meas_stage1)

        # measurement stage 2
        me_circuit.measure(top2reg, meas_stage2)

        return me_circuit

class DiscreteLogMoscaEkertSemiClassicalQFT(DiscreteLogMoscaEkert):
    def __init__(self,
                 b: int,
                 g: int,
                 p: int,
                 r: int = -1,
                 n: int = -1,
                 mod_exp_constructor: Callable[[int, int, int], QuantumCircuit] = None,
                 full_run: bool = False,
                 quantum_instance: Optional[
                     Union[QuantumInstance, Backend, BaseBackend, Backend]] = None) -> None:
        super().__init__(b, g, p, r, n, mod_exp_constructor, full_run, quantum_instance)

    def _exec_qc(self, n, b, g, p, r, mod_exp_constructor, quantum_instance, shots) -> [(int, int, int)]:
        me_circuit = transpile(self._construct_circuit(n, b, g, p, r, mod_exp_constructor),
                               quantum_instance.backend)
        counts = quantum_instance.backend.run(me_circuit, shots=shots).result().get_counts(me_circuit)

        res = list()
        for result in counts.keys():
            # split result measurements (is split bit by bit in this version)
            result_s = result.split(" ")

            #print("Result: ", result)

            # starts with stage 2
            m_stage2_bin = ""
            for i in range(0, n):
                m_stage2_bin = m_stage2_bin + result_s[i]

            m_stage1_bin = ""
            for i in range(0, n):
                m_stage1_bin = m_stage1_bin + result_s[n+i]

            #print(m_stage2_bin)
            #print(m_stage1_bin)

            m_stage1 = int(m_stage1_bin, 2)
            m_stage2 = int(m_stage2_bin, 2)

            res.append((m_stage1, m_stage2, counts[result]))

        return res

    @staticmethod
    def iqft_semi_classical(n: int, circ: QuantumCircuit, reg: QuantumRegister, cllist: [ClassicalRegister]):
        for i in reversed(range(0, n)):  # range is reversed as the first step of IQFT is a swap of the order
            swapped_i = n - i - 1  # index of this qubit after swap would have been performed

            # inverse rotations controlled by less significant qubits (classical now)
            # the rotations are performed with j's in the classical register - those are already implicitly swapped
            # therefore the swapped index is used here
            for j in range(0, swapped_i):
                # The register this qubit was measured in
                clreg = cllist[j]
                circ.p(-np.pi / (2 ** (swapped_i - j)), reg[i]).c_if(clreg, 1)

            circ.h(reg[i])

            # measure it in corresponding classical reg
            circ.measure(reg[i], cllist[swapped_i])

    def _construct_circuit(self, n: int, b: int, g: int, p: int, r: int, mod_exp_constructor) -> QuantumCircuit:
        # infer size of circuit from modular exponentiation circuit
        mod_exp_g = mod_exp_constructor(n, g, p)
        mod_exp_b = mod_exp_constructor(n, b, p)

        iqft = QFT(n).inverse()

        total_circuit_qubits = mod_exp_g.num_qubits
        bottom_register_qubits = total_circuit_qubits - n

        topreg = QuantumRegister(n, "top")
        botreg = QuantumRegister(bottom_register_qubits, "bot")

        cl_stage1 = []  # they have to be adressed as individual classical registers in c_if
        cl_stage2 = []

        me_circuit = QuantumCircuit(topreg, botreg)

        # init measurement registers as collection of single qubit registers
        for i in range(0, n):
            cl1 = ClassicalRegister(1, "f%d" % i)
            cl2 = ClassicalRegister(1, "s%d" % i)

            cl_stage1.append(cl1)
            cl_stage2.append(cl2)

        # Add measurement qubits per stage
        for i in range(0, n):
            me_circuit.add_register(cl_stage1[i])
        for i in range(0, n):
            me_circuit.add_register(cl_stage2[i])

        # First stage

        # H on top
        me_circuit.h(topreg)

        # 1 on bottom
        me_circuit.x(botreg[0])

        # mod exp g^x mod p
        me_circuit.append(mod_exp_g, list(topreg) + list(botreg))

        # semiclassical iqft top register into cl_stage1
        self.iqft_semi_classical(n, me_circuit, topreg, cl_stage1)

        # reset top
        me_circuit.reset(topreg)

        # second stage

        # h on top2
        me_circuit.h(topreg)

        # mod exp b^x' mod p
        me_circuit.append(mod_exp_b, list(topreg) + list(botreg))

        # semiclassical iqft top register into cl_stage2
        self.iqft_semi_classical(n, me_circuit, topreg, cl_stage2)

        return me_circuit

class DiscreteLogMoscaEkertOneQubitQFT(DiscreteLogMoscaEkert):
    def __init__(self,
                 b: int,
                 g: int,
                 p: int,
                 r: int = -1,
                 n: int = -1,
                 mod_mul_constructor: Callable[[int, int, int], QuantumCircuit] = None,
                 full_run: bool = False,
                 quantum_instance: Optional[
                     Union[QuantumInstance, Backend, BaseBackend, Backend]] = None) -> None:
        """
        Needs direct access to the function construction a modular multiplication
        """
        super().__init__(b, g, p, r, n, mod_mul_constructor, full_run, quantum_instance)

        if mod_mul_constructor is None:
            self._mod_exp_constructor = mod_mul_brg
        else:
            self._mod_exp_constructor = mod_mul_constructor

    def _exec_qc(self, n, b, g, p, r, mod_mul_constructor, quantum_instance, shots) -> [(int, int, int)]:
        me_circuit = transpile(self._construct_circuit(n, b, g, p, r, mod_mul_constructor),
                               quantum_instance.backend)
        counts = quantum_instance.backend.run(me_circuit, shots=shots).result().get_counts(me_circuit)

        res = list()
        for result in counts.keys():
            # split result measurements (is split bit by bit in this version)
            result_s = result.split(" ")

            #print("Result: ", result)

            # starts with stage 2
            m_stage2_bin = ""
            for i in range(0, n):
                m_stage2_bin = m_stage2_bin + result_s[i]

            m_stage1_bin = ""
            for i in range(0, n):
                m_stage1_bin = m_stage1_bin + result_s[n+i]

            #print(m_stage2_bin)
            #print(m_stage1_bin)

            m_stage1 = int(m_stage1_bin, 2)
            m_stage2 = int(m_stage2_bin, 2)

            res.append((m_stage1, m_stage2, counts[result]))

        return res

    @staticmethod
    def one_qubit_qft_stage(n: int, a: int, p: int,
                            me_circuit: QuantumCircuit, topreg_single_qubit: QuantumRegister, botreg: QuantumRegister,
                            cllist: [ClassicalRegister], mod_mul_constructor):
        # act like there is a whole register while there actually is only this one qubit
        # this makes the code more similar to the other versions
        topreg = [topreg_single_qubit[0] for i in range(0, n)]

        for i in reversed(range(0, n)):
            # init control to H
            me_circuit.h(topreg[i])

            # exec gate
            mulgate = mod_mul_constructor(n, a**(2**i) % p, p) # U_g^(2^i)

            me_circuit.append(mulgate, [topreg[i], *botreg])

            swapped_i = n - i - 1  # index of this qubit after swap would have been performed
            # inverse rotations controlled by less significant qubits (classical now)

            # the rotations are performed with j's in the classical register - those are already implicitly swapped
            # therefore the swapped index is used here
            for j in range(0, swapped_i):
                # The register this qubit was measured in
                clreg = cllist[j]
                me_circuit.p(-np.pi / (2 ** (swapped_i - j)), topreg[i]).c_if(clreg, 1)

            me_circuit.h(topreg[i])

            # measure it in corresponding classical reg
            me_circuit.measure(topreg[i], cllist[swapped_i])

            # reset to 0 so the rest of the algorithm can reuse this qubit
            me_circuit.reset(topreg[i])

            me_circuit.barrier()  # for visualisation purposes

    def _construct_circuit(self, n: int, b: int, g: int, p: int, r: int, mod_mul_constructor) -> QuantumCircuit:
        # infer size of circuit from modular exponentiation circuit
        mod_mul_dummy = mod_mul_constructor(n, g, p)

        iqft = QFT(n).inverse()

        total_circuit_qubits = mod_mul_dummy.num_qubits  #mod mul has botreg + 1 qubits already
        bottom_register_qubits = total_circuit_qubits - 1

        topreg = QuantumRegister(1, "top")
        botreg = QuantumRegister(bottom_register_qubits, "bot")

        cl_stage1 = []  # they have to be adressed as individual classical registers in c_if
        cl_stage2 = []

        me_circuit = QuantumCircuit(topreg, botreg)

        # init measurement registers as collection of single qubit registers
        for i in range(0, n):
            cl1 = ClassicalRegister(1, "f%d" % i)
            cl2 = ClassicalRegister(1, "s%d" % i)

            cl_stage1.append(cl1)
            cl_stage2.append(cl2)

        # Add measurement qubits per stage
        for i in range(0, n):
            me_circuit.add_register(cl_stage1[i])
        for i in range(0, n):
            me_circuit.add_register(cl_stage2[i])

        # First stage
        me_circuit.x(botreg[0])

        # performs stage1 = apply gate to superposition and do iqft in one qubit
        self.one_qubit_qft_stage(n, g, p, me_circuit, topreg, botreg, cl_stage1, mod_mul_constructor)

        # Second stage
        self.one_qubit_qft_stage(n, b, p, me_circuit, topreg, botreg, cl_stage2, mod_mul_constructor)

        return me_circuit