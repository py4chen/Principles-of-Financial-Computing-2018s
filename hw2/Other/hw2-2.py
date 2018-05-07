#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from math import *
import logging
import numpy as np

logger = logging.getLogger('logger')
logger.setLevel(logging.WARNING)
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
logger.addHandler(ch)

def BOPF():

    # Initial value
    # S = 100, X = 100, H = 110, T = 1, r = 0.05, s = 0.30, and n = 200, and k = 100
    S = 100
    X = 100
    t = 1
    n = 200
    H = 110
    s = 0.3
    r = 0.05
    k = 100
    u = exp(s * sqrt(t / n))
    d = 1 / u   # d = exp(-s * sqrt(t / n))
    r_ = r * t / n
    R = exp(r_)
    p = (R - d) / (u - d)  # Risk-neutral P

    def Amax(j, i):
        maxsum = (S * ((1 - u ** (j - i + 1)) / (1 - u) + u ** (j-i) * d * (1 - d ** i) / (1 - d) ))
        return maxsum / (j + 1)

    def Amin(j, i):
        minsum = (S * ((1 - d ** (i + 1)) / (1 - d) + d ** i * u * (1 - u ** (j - i)) / (1 - u) ))
        return minsum / (j + 1)

    def Average(m, j, i):
        return (((k - m) / k) * Amin(j, i) + (m / k) * Amax(j, i))

    def findl(A, j, i):
        logger.debug('l not found')
        if A < Average(0, j, i):
            # logger.warning('l return 0')
            return 0
        if (A >= Average(k, j, i)):
            # logger.warning('l return k')
            return k
        for l in range(k):
            if Average(l, j, i) <= A and A <= Average(l+1, j, i):
                return l
        # logger.warning('l not found')
        return 0

    # Vanilla
    V = [[max(0, Average(m, n, i) - X) * (Average(m, n, i) < H) for m in range(k + 1)] for i in range(n + 1)]
    V = np.array(V)

    # Barrier
    C = [[max(0, Average(m, n, i) - X) * (Average(m, n, i) < H) for m in range(k+1)] for i in range(n+1)]
    C = np.array(C)
    logger.debug("C: %r", C)

    D = [None] * (k+1)
    D = np.array(D)


    for j in reversed(range(n)):
        logger.warning("loop j=%s", j)
        for i in range(j+1):
            for m in range(k+1):

                # Asian barrier option
                a = Average(m, j, i)
                A_u = ((j+1) * a + S * u ** (j+1-i) * d ** i) / (j+2)
                # logger.debug("A_u: %s", A_u)
                # logger.debug("Amax: %s", Amax(j, i))
                # logger.debug("Amin: %s", Amin(j, i))
                l = findl(A_u, j+1, i)
                try:
                    if l not in [0, k]:
                        x = (A_u - Average(l+1, j+1, i)) / (Average(l, j+1, i) - Average(l+1, j+1, i))
                        C_u = x * C[i][l] + (1-x) * C[i][l+1]
                    else:
                        C_u = C[i][l]
                except:
                    C_u = C[i][l]

                A_d = ((j+1) * a + S * u ** (j-i) * d ** (i+1)) / (j+2)
                l = findl(A_d, j+1, i+1)
                try:
                    if l not in [0, k]:
                        x = (A_d - Average(l+1, j+1, i+1)) / (Average(l, j+1, i+1) - Average(l+1, j+1, i+1))
                        C_d = x * C[i+1][l] + (1-x) * C[i+1][l+1]
                    else:
                        C_d = C[i+1][l]
                except:
                    C_d = C[i+1][l]

                D[m] = 0 if a >= H else ((p * C_u + (1-p) * C_d) / R)

            C[i][:] = D[:]

        for i in range(j + 1):
            for m in range(k + 1):

                # Asian barrier option
                a = Average(m, j, i)
                A_u = ((j + 1) * a + S * u ** (j + 1 - i) * d ** i) / (j + 2)
                # logger.debug("A_u: %s", A_u)
                # logger.debug("Amax: %s", Amax(j, i))
                # logger.debug("Amin: %s", Amin(j, i))
                l = findl(A_u, j + 1, i)
                try:
                    if l not in [0, k]:
                        x = (A_u - Average(l + 1, j + 1, i)) / (Average(l, j + 1, i) - Average(l + 1, j + 1, i))
                        C_u = x * C[i][l] + (1 - x) * C[i][l + 1]
                    else:
                        C_u = C[i][l]
                except:
                    C_u = C[i][l]

                A_d = ((j + 1) * a + S * u ** (j - i) * d ** (i + 1)) / (j + 2)
                l = findl(A_d, j + 1, i + 1)
                try:
                    if l not in [0, k]:
                        x = (A_d - Average(l + 1, j + 1, i + 1)) / (
                        Average(l, j + 1, i + 1) - Average(l + 1, j + 1, i + 1))
                        C_d = x * C[i + 1][l] + (1 - x) * C[i + 1][l + 1]
                    else:
                        C_d = C[i + 1][l]
                except:
                    C_d = C[i + 1][l]

                D[m] = ((p * C_u + (1 - p) * C_d) / R)

            V[i][:] = D[:]

    print("Barrier:", C[0][0])
    print("Vanilla:", V[0][0])
    print("Difference:", V[0][0] - C[0][0])

"""
Barrier: 0.28480283636553444
Vanilla: 0.294358743569828
Difference: 0.009555907204293557
"""

if __name__ == "__main__":
    BOPF()