import numpy as np


A = np.array([1,2,3,4,5])
C = np.array([4,10,18, 28, 40])

# var = np.vstack([A, B, np.ones(len(A))]).T
# var = np.vstack([A, B, np.ones(len(A))]).T
# print(np.linalg.lstsq(var, C, rcond=-1)[0])

print(np.polyfit(A,C, 2))

print(np.random.normal(0, 1, 6))
print(A[-5])
print(np.std(np.array([1,1,1,1])))

D = np.array([[1,2,3],[4,5,6]])
print(np.mean(D[:,:2], axis=1))



rnormal = np.random.normal(0, 1, 1000)
print('[' + ','.join('{}'.format(v) for k, v in enumerate(rnormal)) + ']')