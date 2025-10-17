from math import sqrt
from match_to_py import MSSM

d1 = {
    "g1": 0.15 / sqrt(2),
    "g3": 1.1,
    "cgamma": 0.01,
    "yu11": 0.00001,
    "yu22": 0.007,
    "yu33": 0.9,
    "mPhi": 1001000,
    "muTilde": 1000000,
    "m1": 1.2,
    "m2": 1002000,
    "m3": 1003000,
    "met1": 1004000,
    "met2": 1005000,
    "met3": 1006000,
    "mlt1": 1007000,
    "mlt2": 1008000,
    "mlt3": 1009000,
    "mqt1": 1010000,
    "mqt2": 1011000,
    "mqt3": 1012000,
    "mut1": 1013000,
    "mut2": 1014000,
    "mut3": 2.0,
    "mdt1": 1016000,
    "mdt2": 1017000,
    "mdt3": 1018000,
}

mubarsq = 1.05 * 1.05
hbar = 0.006332574
m1 = MSSM(d1)

print(m1.cH(mubarsq, hbar))
print(m1.cllHH(0, 0, mubarsq, hbar))
print(m1.cHe(1, 1, mubarsq, hbar))
print(m1.cuu(2, 2, 2, 2, mubarsq, hbar))
