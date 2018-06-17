import numpy as np
import copy
import sys
import math


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

        self.ETA = dict()
        self.H2 = dict()
        self.P = dict()
        self.TREE = [list() for _ in range(E + 1)]
        self.TREE_VAR = dict()

    def buildTree(self):
        self.TREE[0].append(0)
        self.TREE_VAR[(0, 0)] = list()
        for k in range(self.n2):
            self.H2[(0, 0, k)] = self.gamma ** 2
            self.TREE_VAR[(0, 0)].append(k)

        # Forward
        for i in range(0, self.E):  # Days
            for j in self.TREE[i]:
                for k in range(self.n2):
                    self.ETA[(i, j, k)] = 0
                for k in range(self.n2):
                    for l in range(2 * self.n1 + 1):
                        self.P[(i, j, k, l)] = 0
                for k in range(self.n2):
                    h2 = self.H2[(i, j, k)]
                    # Find eta
                    eta = int(np.ceil((h2 ** 0.5) / self.gamma))
                    for e in range(eta, sys.maxsize):
                        low = abs(self.r - (h2 / 2.)) / (2. * eta * self.gamma * (self.n1 ** .5))
                        high = min(1 - low, .5)
                        mid = h2 / (2. * eta * eta * self.gamma * self.gamma)
                        if low <= mid <= high:
                            eta = e
                            break
                    if eta == 0:
                        continue

                    self.ETA[(i, j, k)] = eta
                    pu = h2 / (eta * eta * self.gamma * self.gamma) * 0.5 + \
                         (self.r - h2 / 2.) / (2. * eta * self.gamma * (self.n1 ** .5))
                    pm = 1. - h2 / (eta * eta * self.gamma * self.gamma)
                    pd = h2 / (eta * eta * self.gamma * self.gamma) * 0.5 - \
                        (self.r - h2 / 2.) / (2. * eta * self.gamma * (self.n1 ** .5))
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
                    cols = [0] * (2 * n1 + 1)
                    for f in factors:
                        tmp = 1
                        counter = 0
                        for f2 in f:
                            tmp *= umd[f2]
                            counter += f2
                        cols[counter] += tmp

                    for l in range(-self.n1, self.n1 + 1):
                        self.P[(i, j, k, l + self.n1)] = cols[self.n1 - l]

            # max & min of variances
            for j in self.TREE[i]:
                for k in range(self.n2):
                    eta = self.ETA[(i, j, k)]
                    if eta == 0:
                        continue
                    h2 = self.H2[(i, j, k)]

                    for l in range(-self.n1, self.n1 + 1):
                        j_next = j + eta * l
                        epsilon = (l * eta * self.gamma_n - self.r + h2 / 2.) / (h2 ** 0.5)
                        h2_next = self.b0 + self.b1 * h2 + self.b2 * h2 * (epsilon - self.c) ** 2.
                        if j_next not in self.TREE[i + 1]:
                            self.TREE[i + 1].append(j_next)
                            self.TREE_VAR[(i + 1, j_next)] = list()
                            for k in range(self.n2):
                                self.H2[(i + 1, j_next, k)] = h2_next
                                self.TREE_VAR[(i + 1, j_next)].append(k)

                        min_index = min(self.TREE_VAR[(i + 1, j_next)])
                        max_index = max(self.TREE_VAR[(i + 1, j_next)])
                        min_ = self.H2[(i + 1, j_next, min_index)]
                        max_ = self.H2[(i + 1, j_next, max_index)]
                        self.H2[(i + 1, j_next, min_index)] = min(h2_next, min_)
                        self.H2[(i + 1, j_next, max_index)] = max(h2_next, max_)

            # interpolation
            for j_next in self.TREE[i+1]:
                min_index = min(self.TREE_VAR[(i + 1, j_next)])
                max_index = max(self.TREE_VAR[(i + 1, j_next)])
                h2_min = self.H2[(i + 1, j_next, min_index)]
                h2_max = self.H2[(i + 1, j_next, max_index)]
                for mid in range(1, self.n2 - 1):
                    self.H2[(i + 1, j_next, mid)] = h2_min + mid * (h2_max - h2_min) / float(self.n2 - 1)

    def price(self):
        PRICE = dict()
        for j in self.TREE[-1]:
            put = max(self.X - self.S * np.exp(self.gamma_n * j), 0.)
            for k in range(self.n2):
                PRICE[(self.E, j, k)] = put

        # backward
        for i in range(self.E - 1, -1, -1):
            for j in self.TREE[i]:
                for k in range(self.n2):
                    eta = self.ETA[(i, j, k)]
                    if eta == 0:
                        continue
                    h2 = self.H2[(i, j, k)]

                    tmp = 0.
                    for l in range(-self.n1, self.n1 + 1):
                        j_next = j + eta * l
                        eps = (l * eta * self.gamma_n - self.r + h2 / 2.) / (h2 ** .5)
                        h2_next = self.b0 + self.b1 * h2 + self.b2 * h2 * ((eps - self.c) ** 2.)
                        for k_next in range(self.n2 - 1):
                            low = self.H2[(i + 1, j_next, k_next)]
                            high = self.H2[(i + 1, j_next, k_next + 1)]
                            if low <= h2_next <= high:
                                break
                        x = (high - h2_next) / (high - low) if high - low != 0 else 0
                        put_next = x * PRICE[(i + 1, j_next, k_next)] + \
                               (1. - x) * PRICE[(i + 1, j_next, k_next + 1)]
                        tmp += self.P[(i, j, k, l + self.n1)] * put_next

                    PRICE[(i, j, k)] = max(tmp / np.exp(self.r), 0)

        return PRICE[(0, 0, 0)]


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

    r = r * .01  # %

    model = EuropeanPutGarchModel(E, r, S, h0, b0, b1, b2, c, X, n1, n2)
    model.buildTree()
    price = model.price()
    print("PRICE: ", price)


