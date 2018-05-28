import math
import numpy as np
np.random.seed(seed=12345)


def calAmericanStyleAsianPuts(Spot, sigma, n, m, Strike, r, nT):

    T = nT / n
    priceMatrix = np.zeros((m, n))
    for i in range(m):
        priceMatrix[i, :] = Spot * np.exp(
            np.cumsum(
                (r - 0.5 * sigma ** 2) * T
                + (sigma * math.sqrt(T) * np.random.normal(0, 1, n))
            )
        )

    # Asian
    meanMatrix = np.concatenate((np.reshape(priceMatrix[:, 0], (m, -1)), np.zeros((m, n - 1))), axis=1)
    for i in range(1, n):
        meanMatrix[:, i] = np.mean(priceMatrix[:, :i + 1], axis=1)

    diffMatrix = np.maximum(0, Strike - meanMatrix)
    X = np.where(diffMatrix > 0, priceMatrix, 0)
    Xsh = X[:, :-1]
    Y1 = diffMatrix * math.exp(-r * T)
    Y2 = np.concatenate((np.zeros((m, n - 1)), np.vstack(Y1[:, n - 1])), axis=1)
    CV = np.zeros((m, n - 1))

    # iteration
    for i in range(n - 2, -1, -1):
        degree = 5
        reg1 = np.polyfit(Xsh[:, i], Y2[:, i + 1], degree)
        for d in range(degree + 1):
            CV[:, i] += reg1[-d - 1] * Xsh[:, i] ** d
            # CV[:, i] = reg1[2] + reg1[1] * Xsh[:,i] + reg1[0] * (Xsh[:,i] ** 2)

        CV[:, i] = np.nan_to_num(CV[:, i])
        Y2[:, i] = np.where(diffMatrix[:, i] > CV[:, i], Y1[:, i], Y2[:, i + 1] * math.exp(-r * T))

    CV = np.nan_to_num(CV)
    CVp = np.concatenate((CV, np.zeros((m, 1))), axis=1)
    pofMatrix = np.where(CVp > diffMatrix, 0, diffMatrix)

    # first value row
    M = np.zeros((m,n))
    for i in range(m):
        M[i, :] = np.cumsum(pofMatrix[i, :])
    M2 = np.concatenate((np.zeros((m, 1)), M[:, :-1]), axis=1)
    fpofMatrix = np.zeros((m, n))
    for i in range(pofMatrix.shape[1]):
        fpofMatrix[:, i] = np.where(M2[:, i] > 0, 0, pofMatrix[:, i])

    dfpofMatrix = np.zeros((m, n))
    for i in range(n):
        dfpofMatrix[:, i] = fpofMatrix[:, i] * math.exp(-T * r * (i + 1))
    PRICE = np.mean(np.sum(dfpofMatrix, axis=1))
    STD = np.std(np.sum(dfpofMatrix, axis=1)) / (m ** 0.5)

    return PRICE, STD


if __name__ == "__main__":

    isAuto = input("Calculate example? (Y/N)")
    if(len(isAuto) == 0 or isAuto[0] == "Y" or isAuto[0] == "y"):
        # Test Case
        S = 100
        X = 100
        T = 1.
        r = 0.05
        s = 0.30
        n = 100
        k = 100000

        print("S:", S)  # spot price
        print("X:", X)  # strike price
        print("T:", T)  # years
        print("r:", r)  # risk-free interest rate
        print("s:", s)  # volatility
        print("n:", n)  # number of periods
        print("k:", k)  # number of simulation paths


    else:
        S = float(input("S:"))
        X = float(input("X:"))
        T = float(input("T:"))
        r = float(input("r:"))
        s = float(input("s:"))
        n = int(input("n:"))
        k = int(input("k:"))

    print("~~~~~~Start Calculate:~~~~~~")
    price, std = calAmericanStyleAsianPuts(Spot=S, sigma=s, n=n, m=k, Strike=X, r=r, nT=T)
    print("Price:", price)
    print("standard error:", std)