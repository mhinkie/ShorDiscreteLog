# ShorDiscreteLog
Implementation of Shor's Algorithm for the DLP.

The implementation of the algorithm itself is based on
"The Hidden Subgroup Problem and Eigenvalue Estimation on a Quantum Computer" (Michele Mosca, Artur Ekert)
and the Appendix of 
"A Quantum “Magic Box” for the Discrete Logarithm Problem" (Burton S. Kaliski Jr.)

## Modular Exponentiation
The repository contains three gates performing modular exponentiation:

### shor.arithmetic.brg_mod_exp.mod_exp_brg 
Implementation using QFT adder. See

"Circuit for Shor’s algorithm using 2n+3 qubits" (Stéphane Beauregard) 

and 

"Addition on a Quantum Computer" (Thomas G. Draper)

### shor.arithmetic.hrs_mod_exp.mod_exp_hrs
Implementation using bitwise carry addition. See

"Factoring using 2n + 2 qubits with Toffoli based modular multiplication" (Thomas Häner, Martin Roetteler, Krysta M. Svore)

### shor.arithmetic.montgomery_mod_exp.mod_exp_montgomery
Implementation using QFT adder and montgomery multiplication. See

"High Performance Quantum Modular Multipliers" (Rich Rines, Isaac Chuang)

## The Algorithm

There are two versions of the algorithm, each able to use all of the modular exponentiation circuits mentioned above.

The first version uses two seperate registers for the inputs of the two exponentiation operations and measures once at the end of the circuit:

```python
from qiskit import Aer
from shor.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertSeperateRegister

simulator = Aer.get_backend('aer_simulator')

me_algo = DiscreteLogMoscaEkertSeperateRegister(b=4, g=2, p=17, quantum_instance=simulator)
result = me_algo.run()

print(result)
```

This will find `m` such that `g**m % p == 4`. It will also return a success probability in the result, although this value only contains meaningful content 
if `full_run=True` is set in the constructor. When executing a full run, the algorithm will not stop if the right result is found, but will check how many 
shots of the quantum algorithm returned this right result and how many failed.

The second version uses only one register as the input for each stage. This register will be reset after the first measurement:

```python
from qiskit import Aer
from shor.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertSharedRegister

simulator = Aer.get_backend('aer_simulator')

me_algo = DiscreteLogMoscaEkertSharedRegister(b=4, g=2, p=17, full_run=True, quantum_instance=simulator)
result = me_algo.run()

print(result)
```

Per default, the algorithm will use `shor.arithmetic.brg_mod_exp.mod_exp_brg` as the implementation of modular exponentation. It does however accept a function 
which returns any other implementation as a QuantumCircuit and will use this circuit:

```python
from qiskit import Aer, QuantumCircuit

from shor.arithmetic.hrs_mod_exp import mod_exp_hrs
from shor.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertSharedRegister

simulator = Aer.get_backend('aer_simulator')

def mod_exp_qc(n: int, a: int, p: int) -> QuantumCircuit:
    return mod_exp_hrs(n, n, a, p)

me_algo = DiscreteLogMoscaEkertSharedRegister(b=4, g=2, p=17, mod_exp_constructor=mod_exp_qc, quantum_instance=simulator)
result = me_algo.run()

print(result)
```

In `mod_exp_qc`, `n` is expected to be the size of the input register and the circuit is expected to calculate `a**x % p`.

Refer to 
 - [mod_exp_brg](shor/mosca_ekert/me_brg.py)
 - [mod_exp_hrs](shor/mosca_ekert/me_hrs.py)
 - [mod_exp_montgomery](shor/mosca_ekert/me_montgomery.py)

on how the different implementations can be used.

### QFT optimizations
There are two versions of the algorithm, that optimize the the inverse QFT required at the end of each stage.

The first one performs the inverse QFT semi-classically, by conditioning the rotations not on quantum register, but measured values in classical registers:
```python
from qiskit import Aer, QuantumCircuit

from shor.arithmetic.hrs_mod_exp import mod_exp_hrs
from shor.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertSemiClassicalQFT

simulator = Aer.get_backend('aer_simulator')

def mod_exp_qc(n: int, a: int, p: int) -> QuantumCircuit:
    return mod_exp_hrs(n, n, a, p)

me_algo = DiscreteLogMoscaEkertSemiClassicalQFT(b=4, g=2, p=17, quantum_instance=simulator)
result = me_algo.run()

print(result)
```

The second one restructures the algorithm such that only one control register is required (see Section 5.2 in the references paper by Mosca and Ekert).
```python
from qiskit import Aer, QuantumCircuit

from shor.arithmetic.hrs_mod_exp import mod_exp_hrs
from shor.mosca_ekert.mosca_ekert import DiscreteLogMoscaEkertOneQubitQFT

simulator = Aer.get_backend('aer_simulator')

me_algo = DiscreteLogMoscaEkertOneQubitQFT(b=4, g=2, p=17, quantum_instance=simulator)
result = me_algo.run()

print(result)
```

Note that in this version the algorithm uses multiplication gates directly (not by using exponentiation gates). This means that not the implementation of the 
modular exponentiation gate can be manually supplied but instead the implementation of the modular multiplication gate as `mod_mul_constructor`.

#### Requirements
Qiskit, Numpy for the implementations. Sympy and Pytest for the tests.
