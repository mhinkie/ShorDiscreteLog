from qiskit import Aer

from shor.arithmetic.montgomery_mod_exp import mod_exp_montgomery
from shor.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertSharedRegister

simulator = Aer.get_backend('aer_simulator')

def get_order(g, p):
    for x in range(1, p):
        if g**x % p == 1:
            return x

failed_attempts = list()
lowest_attempt_prob = 100
lowest_attempt_vals = None

p = 17
for g in range(2, p):
    # find period r
    r = get_order(g, p)

    if r >= 0:
        for m in range(0, r):
            b = g**m % p
            print("g=%d, r=%d, b=%d, exp_m=%d" % (g, r, b, m))
            me_algo = DiscreteLogMoscaEkertSharedRegister(b, g, p, r=r, full_run=True, quantum_instance=simulator, mod_exp_constructor=mod_exp_montgomery)
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
