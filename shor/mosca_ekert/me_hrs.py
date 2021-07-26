from arithmetic.brg_mod_exp import mod_exp_brg
from qiskit.circuit.library import QFT
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, assemble, Aer, transpile, execute
from qiskit.circuit.library import XGate
import numpy as np
import math

from arithmetic.hrs_mod_exp import mod_exp_hrs
from ma_impl.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertSharedRegister
from util import set_input

simulator = Aer.get_backend('aer_simulator')

def get_order(g, p):
    for x in range(1, p):
        if g**x % p == 1:
            return x

def mod_exp_constr_hrs(n, a, p):
    return mod_exp_hrs(n, n, a, p)

failed_attempts = list()
lowest_attempt_prob = 100
lowest_attempt_vals = None

p = 37
for g in range(2, p):
    # find period r
    r = get_order(g, p)

    if r >= 0:
        for m in range(0, r):
            b = g**m % p
            print("g=%d, r=%d, b=%d, exp_m=%d" % (g, r, b, m))
            me_algo = DiscreteLogMoscaEkertSharedRegister(b, g, p, r=r, full_run=True, quantum_instance=simulator, mod_exp_constructor=mod_exp_constr_hrs)
            res = me_algo.run()

            if res["success_prob"] == 0.00:
                failed_attempts.append((g, r, b, m, res))
            else:
                if res["success_prob"] < lowest_attempt_prob:
                    lowest_attempt_vals = (g, r, b, m)
                    lowest_attempt_prob = res["success_prob"]

print("Failed attempts: ")
print(failed_attempts)
print("Lowest prob: ")
print(lowest_attempt_prob)
print(lowest_attempt_vals)
