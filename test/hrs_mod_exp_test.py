import unittest
from qiskit import QuantumCircuit, QuantumRegister, Aer

from shor.arithmetic.hrs_mod_exp import mod_exp_hrs
from .test_util import run_on_simulator, set_input
from sympy.ntheory.generate import prevprime

import random

def split_result(result_str, n):
    """returns (top_reg, bot_reg, anc)"""
    return result_str[-n:], result_str[-2 * n:-n], result_str[0:-2 * n]

class HRSModExpTest(unittest.TestCase):
    def setUp(self):
        self.simulator = Aer.get_backend('aer_simulator')

    def test_mod_exp(self):
        # Perform modular exponentiation 10 times for 1024 shots each
        # using H on the input |x>, while setting a random value
        # for the content of the second register and using
        # a different randomly selected input constant each time
        # the module is always chosen to be the largest prime < 2^n

        n = 4 # bits for input

        for _ in range(0, 10):
            p = prevprime(2**n)
            input_constant = random.randint(1, p-1)
            bot_reg_content = random.randint(0, p-1)

            print("Input constant:", input_constant)
            print("Bot reg content: ", bot_reg_content)
            print("Module p=", p)

            # set up circuit
            top = QuantumRegister(n, "top")
            bot = QuantumRegister(2 * n + 1, "bot")

            circ = QuantumCircuit(top, bot)

            circ.h(top)
            circ.append(set_input(n, bot_reg_content), bot[:n])

            circ.append(mod_exp_hrs(n, n, input_constant, p), circ.qubits)

            circ.measure_all()

            counts = run_on_simulator(circ, self.simulator, 1024)

            for result in counts.keys():
                (top_result, bot_result, anc_result) = split_result(result, n)

                # ancilla qubits should be 0
                self.assertEqual(int(anc_result, 2), 0, "ancilla content should be zero, got result: %s, for input %s" %
                                 (result, str((p, input_constant, bot_reg_content))))

                # bot reg content should be bot * a^(top) mod p
                expected_result = (bot_reg_content * (input_constant**int(top_result, 2) % p)) % p
                self.assertEqual(int(bot_result, 2), expected_result, "Expected result %d in bottom register - "
                                                                      "got %d, for input %s - full result %s" %
                                 (expected_result, int(bot_result, 2), str((p, input_constant, bot_reg_content)),
                                  result))



