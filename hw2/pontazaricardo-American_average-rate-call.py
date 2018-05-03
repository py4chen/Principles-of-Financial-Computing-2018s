"""
https://github.com/pontazaricardo/Finance_American_average-rate_call
"""

from math import exp, sqrt
import numpy as np

def findA(m,j,i,S,k,u,d):
    aminji = (1 / (j+1)) * (S * ((1 - d ** (i+1)) / (1-d)) +
                            S * (d ** i) * u * ((1 - u ** (j-i)) / (1-u)))
    aminTerm = ((k - m) / k) * aminji
    amaxji = (1 / (j+1) ) * (S * ((1 - (u ** (j-i+1))) / (1 - u)) +
                             S * (u ** (j-i)) * d * ((1 - d ** i) / (1 - d)))
    amaxTerm = (m / k) * amaxji
    amji = aminTerm + amaxTerm
    return amji


def findLUP(Au,j,i,S,k,u,d):
    epsilon = 1e-04
    lvalue = 0
    x1 = 1
    for l in range(0, k):
        AL1 = findA(l, j+1, i, S, k, u, d)
        AL2 = findA(l+1, j+1, i, S, k, u, d)
        if ((AL1 - epsilon <= Au) and (Au <= AL2+epsilon)):
            lvalue=l
            if (AL1 == AL2):
               x1 = 1
            else:
               x1 = (Au-AL2) / (AL1-AL2)
            break
    return lvalue, x1

def findLDOWN(Ad, j, i, S, k, u, d):
    lvalue = 0
    x2 = 1
    epsilon = 1e-04
    for l in range(0, k):
        AL1 = findA(l, j+1, i+1, S, k, u, d)
        AL2 = findA(l+1, j+1, i+1, S, k, u, d)
        # print(AL1 - epsilon, Ad, (AL1 - epsilon <= Ad), Ad, AL2 + epsilon, (Ad <= AL2 + epsilon))
        if ((AL1 - epsilon <= Ad) and (Ad <= AL2 + epsilon)):
            lvalue = l
            if (AL1 == AL2):
                x2 = 1
            else:
                x2 = (Ad - AL2) / (AL1 - AL2)
            break
    return lvalue, x2



if __name__ == "__main__":

    S = 100
    X = 70
    t = 2
    sigma = 0.2
    r1 = 0.05
    n = 40
    k = 5


    deltaT = t / n
    r = r1 * deltaT
    u = exp(sigma * sqrt(deltaT))
    d = 1 / u
    p = (exp(r) - d) / (u - d)

    C = np.zeros((n+1, k+1))
    D = np.zeros(k+1)
    delta = 0

    for i in range(0, n + 1):
        for m in range(0, k + 1):
            C[i, m] = max(0, findA(m, n, i, S, k, u, d) - X)


    for j in range(n-1, -1, -1):
        for i in range(0, j + 1):
            for m in range(0, k + 1):
                a = findA(m, j, i, S, k, u, d)  # Average

                # Up calculations
                Au = ((j+1) * a + S * (u ** (j+ 1- i)) * (d ** i) ) / (j + 2)
                lup, x1 = findLUP(Au, j, i, S, k, u, d)
                Cu = x1 * C[i, lup] + (1-x1) * C[i, lup+1]   # Linear interpolation

                # Down calculations
                Ad = ((j+1) * a + S * (u ** (j-i)) * (d ** (i+1)))/(j+2)
                ldown, x2 = findLDOWN(Ad, j, i, S, k, u, d)
                Cd = x2 * C[i + 1, ldown] + (1 - x2) * C[i+1, ldown+1]   # Linear interpolation

                # Calculate D
                D[m] = max(a-X, (p * Cu + (1 - p) * Cd) * exp(-r))

                delta = (Cu-Cd) / (S*u-S*d)


            # for w in range(0, k + 1):
            #     C[i,w] = D[w]
            C[i, :] = D

    print('-------------------------')
    print('The price is: ', C[0, 0])
    print('Delta is: ', delta)
    print('-------------------------')


"""
The price is:  36.30801750098361
Delta is:  0.9514529663683023
"""