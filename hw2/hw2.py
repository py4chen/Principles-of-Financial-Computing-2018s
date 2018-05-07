from math import exp, sqrt
import numpy as np
import sys

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

    def calAverage(self, m, j, i):
        if np.isnan(self.A[j, i, m]):
            a_min_ji = (1 / (j + 1)) * (self.S * ((1 - self.d ** (i + 1)) / (1 - self.d)) +
                                      self.S * (self.d ** i) * self.u * ((1 - self.u ** (j - i)) / (1 - self.u)))
            a_min = ((self.k - m) / self.k) * a_min_ji
            a_max_ji = (1 / (j + 1)) * (self.S * ((1 - (self.u ** (j - i + 1))) / (1 - self.u)) +
                                      self.S * (self.u ** (j - i)) * self.d * ((1 - self.d ** i) / (1 - self.d)))
            a_max = (m / self.k) * a_max_ji
            amji = a_min + a_max
            self.A[j, i, m] = amji
        return self.A[j, i, m]

    def calLup(self, Au, j, i):
        epsilon = 1e-04
        lvalue = 0
        x1 = 1
        for l in range(0, self.k):
            AL1 = self.calAverage(l, j + 1, i)
            AL2 = self.calAverage(l + 1, j + 1, i)
            if AL1 - epsilon <= Au <= AL2 + epsilon:
                lvalue = l
                if AL1 == AL2:
                    x1 = 1
                else:
                    x1 = (Au - AL2) / (AL1 - AL2)
                break
        return lvalue, x1

    def calLdown(self, Ad, j, i):
        lvalue = 0
        x2 = 1
        epsilon = 1e-04
        for l in range(0, self.k):
            AL1 = self.calAverage(l, j + 1, i + 1)
            AL2 = self.calAverage(l + 1, j + 1, i + 1)
            # print(AL1 - epsilon, Ad, (AL1 - epsilon <= Ad), Ad, AL2 + epsilon, (Ad <= AL2 + epsilon))
            if AL1 - epsilon <= Ad <= AL2 + epsilon:
                lvalue = l
                if AL1 == AL2:
                    x2 = 1
                else:
                    x2 = (Ad - AL2) / (AL1 - AL2)
                break
        return lvalue, x2

    def price(self):

        # Process
        for i in range(0, self.n + 1):
            for m in range(0, self.k + 1):
                self.C[i, m] = max(0, self.calAverage(m, self.n, i) - self.X)
                self.B[i, m] = max(0, self.calAverage(m, self.n, i) - self.X)

        D = np.zeros(self.k + 1)
        Db = np.zeros(self.k + 1)

        for j in range(self.n - 1, -1, -1):
            sys.stdout.write('\r j = %d      ' % j)
            sys.stdout.flush()

            for i in range(0, j + 1):
                for m in range(0, self.k + 1):
                    a = self.calAverage(m, j, i)  # Average

                    # Up calculations
                    Au = ((j + 1) * a + self.S * (self.u ** (j + 1 - i)) * (self.d ** i)) / (j + 2)
                    lup, x1 = self.calLup(Au, j, i)
                    Cu = x1 * self.C[i, lup] + (1 - x1) * self.C[i, lup + 1]  # Vanilla Linear interpolation
                    Bu = x1 * self.B[i, lup] + (1 - x1) * self.B[i, lup + 1]  # Barrier Linear interpolation

                    # Down calculations
                    Ad = ((j + 1) * a + self.S * (self.u ** (j - i)) * (self.d ** (i + 1))) / (j + 2)
                    ldown, x2 = self.calLdown(Ad, j, i)
                    Cd = x2 * self.C[i + 1, ldown] + (1 - x2) * self.C[i + 1, ldown + 1]  # Vanilla Linear interpolation
                    Bd = x2 * self.B[i + 1, ldown] + (1 - x2) * self.B[i + 1, ldown + 1]  # Barrier Linear interpolation

                    # Calculate D & Db
                    D[m] = (self.p * Cu + (1 - self.p) * Cd) * exp(-self.r)
                    Db[m] = (self.p * Bu + (1 - self.p) * Bd) * exp(-self.r)

                    delta_vanilla = (Cu - Cd) / (self.S * self.u - self.S * self.d)
                    delta_barrier = (Bu - Bd) / (self.S * self.u - self.S * self.d)

                    if a >= self.H:  # Knock-out if average price reaches barrier price
                        Db[m] = 0

                self.C[i, :] = D
                self.B[i, :] = Db


        print('\n-------------------------')
        print('Vanilla price is: ', self.C[0, 0])
        print('Vanilla Delta is: ', delta_vanilla)
        print('-------------------------')
        print('KnockOut price is: ', self.B[0, 0])
        print('KnockOut Delta is: ', delta_barrier)
        print('-------------------------')

        # Knock-in = Vanilla - Knock-out
        differenceValue = self.C[0, 0] - self.B[0, 0]
        print('KnockIn Value is: ', differenceValue)
        print('KnockIn Delta is: ', delta_vanilla - delta_barrier)
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
Calculate example? (Y/N)y
S: 100
X: 100
H: 110
T: 1
r: 0.05
s: 0.3
n: 200
k: 100
~~~~~~Start Calculate:~~~~~~
 j = 0
-------------------------
Vanilla price is:  8.631127489376166
Vanilla Delta is:  0.5637618041777994
-------------------------
KnockOut price is:  0.304825090845017
KnockOut Delta is:  -0.009722096283885115
-------------------------
KnockIn Value is:  8.326302398531148
KnockIn Delta is:  0.5734839004616845
-------------------------

real    11m5.156s
user    11m3.174s
sys     0m0.692s

TA's Correct Value:
call price is 8.3514
its delta is 0.5726
"""