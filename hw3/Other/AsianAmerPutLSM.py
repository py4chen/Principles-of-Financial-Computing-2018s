import math
import numpy as np

Spot = 100
sigma = 0.3
n = 100000
m = 100
Strike = 100
r = 0.05
mT = 1


T = mT / m

np.random.seed(seed=123)

GBM = np.zeros((n, m))

for i in range(n):

    GBM[i, :] = Spot * np.exp(
        np.cumsum(
            (r - 0.5 * sigma ** 2) * T
            + (sigma * math.sqrt(T) * np.random.normal(0, 1, m))
        )
    )

    # GBM[i, :] = Spot * np.exp(
    #     np.cumsum(
    #         (r - 0.5 * sigma ** 2) * T
    #         + (sigma * math.sqrt(T) * np.random.normal(i/n, 0, m))
    #     )
    # )
    #                                                               * rnormal)))

# print(GBM)
CFLp = np.concatenate((np.reshape(GBM[:, 0], (n, -1)), np.zeros((n, m - 1))), axis=1)

for i in range(1, m):
    CFLp[:, i] = np.mean(GBM[:, :i+1], axis=1)


CFL = np.maximum(0, Strike - CFLp)


X = np.where(CFL > 0, GBM, 0)

Xsh = X[:,:-1]

Y1 = CFL * math.exp(-r * T)

Y2 = np.concatenate((np.zeros((n, m - 1)), np.vstack(Y1[:, m-1])), axis=1)

CV = np.zeros((n, m - 1))


for i in range(m-2, -1, -1):
    degree = 5
    reg1 = np.polyfit(Xsh[:,i], Y2[:, i+1], degree)
    for d in range(degree+1):
        CV[:, i] += reg1[-d-1] * Xsh[:,i] ** d
        # CV[:, i] = reg1[2] + reg1[1] * Xsh[:,i] + reg1[0] * (Xsh[:,i] ** 2)

    CV[:, i] = np.nan_to_num(CV[:, i])
    Y2[:, i] = np.where(CFL[:, i] > CV[:, i], Y1[:, i], Y2[:, i + 1] * math.exp(-r * T))

CV = np.nan_to_num(CV)

CVp = np.concatenate((CV, np.zeros((n,1))), axis=1)
POF = np.where(CVp > CFL, 0, CFL)

def firstValueRow(x):
    cumSumMat = np.zeros((x.shape[0], x.shape[1]))
    for i in range(x.shape[0]):
        cumSumMat[i,:] = np.cumsum(x[i,:])
    cumSumMat2 = np.concatenate((np.zeros((x.shape[0],1)), cumSumMat[:,:-1]), axis=1)
    ResultMat = np.zeros((x.shape[0], x.shape[1]))
    for i in range(x.shape[1]):
        ResultMat[:, i] = np.where(cumSumMat2[:, i] > 0, 0, x[:, i])
    return ResultMat

FPOF = firstValueRow(POF)

dFPOF = np.zeros((n,m))
for i in range(m):
    dFPOF[:, i] = FPOF[:, i] * math.exp(-T * r * (i + 1))

PRICE = np.mean(np.sum(dFPOF, axis=1))
STD = np.std(np.sum(dFPOF, axis=1)) / (n ** 0.5)

print(PRICE)
print(STD)




