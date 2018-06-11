
import numpy as np
import copy
from collections import OrderedDict


class EuropeanPutGarchModel:

    def __init__(self, E, rT, S, h0, b0, b1, b2, c, X, n1, n2):
        self.E = E
        self.r = rT / 365.
        self.S = S
        self.gamma = h0
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.c = c
        self.X = X
        self.n1 = n1
        self.n2 = n2

        self.gamma_n = h0 / (n1 ** .5)
        self.eta_tree = [{} for _ in range(E)]
        # Variance tree
        self.h2_tree = [{} for _ in range(E + 1)]  # [ {} * (days + 1) ]
        self.h2_tree[0][0] = [h0 ** 2 for _ in range(n2)] # h2_tree[day][variance]
        # Probability tree
        self.p_tree = [{} for _ in range(E)]

    def buildTree(self):
        # Forward trees buildup using RT algorithm
        for i in range(0, self.E):  # Days
            # Compute eta_tree[ i] and p_tree[i] using h2_tree[i].
            # print(sorted(h2_tree[i].keys()))

            for j in sorted(self.h2_tree[i].keys()):
                self.eta_tree[i][j] = [0 for _ in range(self.n2)]
                self.p_tree[i][j] = [[0 for _ in range(2 * self.n1 + 1)] for _ in range(self.n2)]
                for k in range(self.n2):
                    h2 = self.h2_tree[i][j][k]

                    # Find eta
                    eta = int(np.ceil((h2 ** 0.5) / self.gamma))
                    while True:
                        low_bound = abs(self.r - (h2 / 2.)) / (2. * eta * self.gamma * (self.n1 ** .5))
                        up_bound = min(1 - low_bound, .5)
                        mid = h2 / (2. * eta * eta * self.gamma * self.gamma)
                        if low_bound <= mid <= up_bound:
                            break
                        eta += 1

                    if eta == 0:
                        continue

                    self.eta_tree[i][j][k] = eta

                    # (2 * n1 + 1) coefficients of (pu * x ** 2 + pm * x + pd) ** n1
                    pu = h2 / (eta * eta * self.gamma * self.gamma) * 0.5 + \
                         (self.r - h2 / 2.) / (2. * eta * self.gamma * (self.n1 ** .5))
                    pm = 1. - h2 / (eta * eta * self.gamma * self.gamma)
                    pd = h2 / (eta * eta * self.gamma * self.gamma) * 0.5 - \
                        (self.r - h2 / 2.) / (2. * eta * self.gamma * (self.n1 ** .5))

                    # Calculate Prob
                    factors = [[]]
                    for _ in range(n1):
                        tmp = []
                        for o in [0, 1, 2]:
                            for f in factors:
                                f2 = copy.copy(f)
                                f2.append(o)
                                tmp.append(f2)
                        factors = tmp
                    umd = (pu, pm, pd)
                    coefs = [0] * (2 * n1 + 1)
                    for f in factors:
                        tmp = 1
                        counter = 0
                        for f2 in f:
                            tmp *= umd[f2]
                            counter += f2
                        coefs[counter] += tmp

                    for l in range(-self.n1, self.n1 + 1):
                        self.p_tree[i][j][k][l + self.n1] = coefs[self.n1 - l]
                        # print(eta_tree[i][j])
                        # print(p_tree[i][j])

            # Compute max./min. variances in h2_tree[i + 1].
            for j in sorted(self.h2_tree[i].keys()):
                for k in range(self.n2):
                    eta = self.eta_tree[i][j][k]
                    if eta == 0:
                        continue
                    h2 = self.h2_tree[i][j][k]

                    # get_next_h2
                    for l in range(-self.n1, self.n1 + 1):
                        next_j = j + eta * l
                        epslon = (l * eta * self.gamma_n - self.r + h2 / 2.) / (h2 ** .5)
                        next_h2 = self.b0 + self.b1 * h2 + self.b2 * h2 * (epslon - self.c) ** 2.
                        if next_j not in self.h2_tree[i + 1]:
                            self.h2_tree[i + 1][next_j] = [next_h2 for _ in range(self.n2)]
                        else:
                            min_ = self.h2_tree[i + 1][next_j][0]
                            self.h2_tree[i + 1][next_j][0] = min(next_h2, min_)
                            max_ = self.h2_tree[i + 1][next_j][-1]
                            self.h2_tree[i + 1][next_j][-1] = max(next_h2, max_)

            # Interpolation of variances (for n2 > 2)
            for next_j in sorted(self.h2_tree[i + 1].keys()):
                min_ = self.h2_tree[i + 1][next_j][0]
                max_ = self.h2_tree[i + 1][next_j][-1]
                for k in range(self.n2):
                    self.h2_tree[i + 1][next_j][k] = min_ + k * (max_ - min_) / (self.n2 - 1.)

    def price(self):
        # Pricing at the last day
        put_tree = [{} for _ in range(self.E + 1)]
        for j in sorted(self.h2_tree[-1].keys()):
            # put = max(stock * np.exp(h0 * j) - strike, 0.)
            put = max(self.X - self.S * np.exp(self.gamma_n * j), 0.)
            put_tree[-1][j] = [put for _ in range(self.n2)]

        # Backward induction
        for i in range(self.E - 1, -1, -1):
            for j in sorted(self.h2_tree[i].keys()):
                put_tree[i][j] = [put for _ in range(self.n2)]
                for k in range(self.n2):
                    eta = self.eta_tree[i][j][k]
                    if eta == 0:
                        continue
                    h2 = self.h2_tree[i][j][k]

                    #
                    put = 0.
                    for l in range(-self.n1, self.n1 + 1):
                        next_j = j + eta * l
                        eps = (l * eta * self.gamma_n - self.r + h2 / 2.) / (h2 ** .5)
                        next_h2 = self.b0 + self.b1 * h2 + self.b2 * h2 * ((eps - self.c) ** 2.)

                        # Find the next (k, k+1) interval bounding next_h2.
                        for next_k in range(self.n2 - 1):
                            low = self.h2_tree[i + 1][next_j][next_k]
                            up = self.h2_tree[i + 1][next_j][next_k + 1]
                            if low <= next_h2 <= up:
                                break

                        x = (next_h2 - up) / (low - up) if low - up != 0 else 0
                        put_ = x * put_tree[i + 1][next_j][next_k] + \
                               (1. - x) * put_tree[i + 1][next_j][next_k + 1]

                        put += self.p_tree[i][j][k][l + self.n1] * put_

                    # American put pricing
                    put_tree[i][j][k] = max(put / np.exp(self.r), 0)

        return put_tree[0][0][0]

if __name__ == "__main__":

    isAuto = input("Calculate example? (Y/N)")
    if(len(isAuto) == 0 or isAuto[0] == "Y" or isAuto[0] == "y"):
        # Test Case
        E = 30
        r = 5
        S = 100
        h0 = 0.010469
        b0 = 0.000006575
        b1 = 0.9
        b2 = 0.04
        c = 0
        X = 100
        n1 = 3
        n2 = 3

        print("E:", E)      # days before expiration
        print("r(%):", r)   # interest rate (%)
        print("S:", S)      # stock price at time 0
        print("h0:", h0)    #
        print("b0:", b0)    #
        print("b1:", b1)    #
        print("b2:", b2)    #
        print("c:", c)      #
        print("X:", X)      # strike price
        print("n1:", n1)    # number of partitions per day
        print("n2:", n2)    # number of variances per node


    else:
        E = int(input("E:"))
        r = float(input("r(%):"))
        S = float(input("S:"))
        h0 = float(input("h0:"))
        b0 = float(input("b0:"))
        b1 = float(input("b1:"))
        b2 = float(input("b2:"))
        c = float(input("c:"))
        X = float(input("X:"))
        n1 = int(input("n1:"))
        n2 = int(input("n2:"))

    r = r * .01 # %

    model = EuropeanPutGarchModel(E ,r ,S ,h0,b0,b1,b2,c ,X ,n1,n2)
    model.buildTree()
    price = model.price()
    print("PRICE: ", price)


