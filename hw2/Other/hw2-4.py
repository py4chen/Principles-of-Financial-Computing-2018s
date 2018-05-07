from math import exp, sqrt, floor, log
import numpy as np
import sys

"""
if a > self.H
"""

class EuropeanArithmeticAverageRateKnockInCall:

    def __init__(self, S, X, H, t, sigma, r1, n, k):

        # Initial Value
        self.S = S
        self.X = X
        self.H = H

        deltaT = t / n
        self.r = r1 * deltaT

        self.u = exp(sigma * sqrt(deltaT))
        self.d = 1 / self.u
        self.p = (exp(self.r) - self.d) / (self.u - self.d)

        self.n = n
        self.k = k

        # Vanilla Value
        self.C = np.zeros((n + 1, k + 1))
        self.delta_vanilla = 0

        # Barrier Knock-out Value
        self.B = np.zeros((n + 1, k + 1))
        self.delta_barrier = 0

        # Average, for cache
        self.A = np.empty((n + 1, n + 1, k + 1)) * np.nan

    def findA(self, m, j, i):
        if np.isnan(self.A[j, i, m]):
            aminji = (1 / (j + 1)) * (self.S * ((1 - self.d ** (i + 1)) / (1 - self.d)) +
                                      self.S * (self.d ** i) * self.u * ((1 - self.u ** (j - i)) / (1 - self.u)))
            aminTerm = ((self.k - m) / self.k) * aminji
            amaxji = (1 / (j + 1)) * (self.S * ((1 - (self.u ** (j - i + 1))) / (1 - self.u)) +
                                      self.S * (self.u ** (j - i)) * self.d * ((1 - self.d ** i) / (1 - self.d)))
            amaxTerm = (m / self.k) * amaxji
            amji = aminTerm + amaxTerm
            self.A[j, i, m] = amji
        return self.A[j, i, m]

    def findLUP(self, Au, j, i):
        epsilon = 1e-04
        lvalue = 0
        x1 = 1
        for l in range(0, self.k):
            AL1 = self.findA(l, j + 1, i)
            AL2 = self.findA(l + 1, j + 1, i)
            if ((AL1 - epsilon <= Au) and (Au <= AL2 + epsilon)):
                lvalue = l
                if (AL1 == AL2):
                    x1 = 1
                else:
                    x1 = (Au - AL2) / (AL1 - AL2)
                break
        return lvalue, x1

    def findLDOWN(self, Ad, j, i):
        lvalue = 0
        x2 = 1
        epsilon = 1e-04
        for l in range(0, self.k):
            AL1 = self.findA(l, j + 1, i + 1)
            AL2 = self.findA(l + 1, j + 1, i + 1)
            # print(AL1 - epsilon, Ad, (AL1 - epsilon <= Ad), Ad, AL2 + epsilon, (Ad <= AL2 + epsilon))
            if ((AL1 - epsilon <= Ad) and (Ad <= AL2 + epsilon)):
                lvalue = l
                if (AL1 == AL2):
                    x2 = 1
                else:
                    x2 = (Ad - AL2) / (AL1 - AL2)
                break
        return lvalue, x2

    def price(self):

        # Process
        for i in range(0, self.n + 1):
            for m in range(0, self.k + 1):
                self.C[i, m] = max(0, self.findA(m, self.n, i) - self.X)
                self.B[i, m] = max(0, self.findA(m, self.n, i) - self.X)

        D = np.zeros(self.k + 1)
        Db = np.zeros(self.k + 1)

        for j in range(self.n - 1, -1, -1):
            sys.stdout.write('\r j = %d' % j)
            sys.stdout.flush()


            for i in range(0, j + 1):
                for m in range(0, self.k + 1):
                    a = self.findA(m, j, i)  # Average

                    # Up calculations
                    Au = ((j + 1) * a + self.S * (self.u ** (j + 1 - i)) * (self.d ** i)) / (j + 2)
                    lup, x1 = self.findLUP(Au, j, i)
                    Cu = x1 * self.C[i, lup] + (1 - x1) * self.C[i, lup + 1]  # Vanilla Linear interpolation
                    Bu = x1 * self.B[i, lup] + (1 - x1) * self.B[i, lup + 1]  # Barrier Linear interpolation

                    # Down calculations
                    Ad = ((j + 1) * a + self.S * (self.u ** (j - i)) * (self.d ** (i + 1))) / (j + 2)
                    ldown, x2 = self.findLDOWN(Ad, j, i)
                    Cd = x2 * self.C[i + 1, ldown] + (1 - x2) * self.C[i + 1, ldown + 1]  # Vanilla Linear interpolation
                    Bd = x2 * self.B[i + 1, ldown] + (1 - x2) * self.B[i + 1, ldown + 1]  # Barrier Linear interpolation

                    # Calculate D & Db
                    D[m] = (self.p * Cu + (1 - self.p) * Cd) * exp(-self.r)
                    Db[m] = (self.p * Bu + (1 - self.p) * Bd) * exp(-self.r)

                    delta_vanilla = (Cu - Cd) / (self.S * self.u - self.S * self.d)
                    delta_barrier = (Bu - Bd) / (self.S * self.u - self.S * self.d)

                    if a > self.H:  # Knock-out if average price reaches barrier price
                        Db[m] = 0

                self.C[i, :] = D
                self.B[i, :] = Db


        print('-------------------------')
        print('Vanilla price is: ', self.C[0, 0])
        print('Vanilla Delta is: ', delta_vanilla)
        print('-------------------------')
        print('Barrier E price is: ', self.B[0, 0])
        print('Vanilla Delta is: ', delta_barrier)
        print('-------------------------')

        # Knock-in = Vanilla - Knock-out
        # I'm not sure whether this method to calculate knock-in delta is correct
        differenceValue = self.C[0, 0] - self.B[0, 0]
        print('differenceValue is: ', differenceValue)
        print('differenceDelta is: ', delta_vanilla - delta_barrier)
        print('-------------------------')

if __name__ == "__main__":

    isAuto = input("Calculate example? (Y/N)")
    if(len(isAuto) == 0 or isAuto[0] == "Y" or isAuto[0] == "y"):
        # Test Case
        S = 100  # (1) S (spot price)
        X = 100  # (2) X (strike price)
        H = 110  # (3) H (barrier price)
        T = 1  # (4) T (years)
        r = 0.05  # (5) r (risk-free interest rate)
        s = 0.3  # (6) s (volatility)
        n = 200  # (7) n (number of periods)
        k = 100  # (8) k (number of buckets)

        print("S:", S)  # spotPrice
        print("X:", X)  # strikePrice
        print("H:", H)  # barrier price
        print("T:", T)  # years
        print("r:", r)  # risk-free interest rate
        print("s:", s)  # volatility
        print("n:", n)  # number of periods
        print("k:", k)  # number of buckets


    else:
        S = float(input("S:"))  # spotPrice
        X = float(input("X:"))  # strikePrice
        H = float(input("H:"))  # barrier price
        T = float(input("T:"))  # years
        r = float(input("r:"))  # risk-free interest rate
        s = float(input("s:"))  # volatility
        n = int(input("n:"))  # number of periods
        k = int(input("k:"))  # number of buckets

    print("~~~~~~Start Calculate:~~~~~~")
    EuropeanArithmeticAverageRateKnockInCall(S=S, X=X, H=H, t=T, sigma=s, r1=r, n=n, k=k).price()


"""
Vanilla price is:  8.631127489376166
Vanilla Delta is:  0.5637618041777994
-------------------------
Barrier E price is:  0.304825090845017
Vanilla Delta is:  -0.009722096283885115
-------------------------
differenceValue is:  8.326302398531148
differenceDelta is:  0.5734839004616845
-------------------------

real    10m22.314s
user    10m21.107s
sys     0m0.633s
"""