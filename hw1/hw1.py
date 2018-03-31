import math

class AmericanPutCRRBinomialTree:
    class Node:
        def __init__(self, price):
            self.price = price
            self.value = 0
            self.upParent = None
            self.downParent = None
            self.downChild = None
            self.upChild = None

    def __init__(self, S, K, r, s, T, n):
        self.S = S
        self.K = K
        self.r= r
        self.s = s
        self.T = T
        self.n = n

        self.u = math.exp(s * math.sqrt(self.T / self.n))  # Page 284
        self.d = 1 / self.u
        self.p = (math.exp(self.r * self.T / self.n) - self.d) / (self.u - self.d)

    def price(self):
        self.root = self.Node(self.S)
        self.levels = [[self.root]]

        for period in range(self.n):
            currentLevel = self.levels[-1]
            childLevel = [AmericanPutCRRBinomialTree.Node(price=currentLevel[0].price * self.u)]
            for currentNode in currentLevel:
                currentNode.upChild = childLevel[-1]
                currentNode.upChild.downParent = currentNode
                downChild = AmericanPutCRRBinomialTree.Node(price=currentNode.price * self.d)
                downChild.upParent = currentNode
                currentNode.downChild = downChild
                childLevel.append(downChild)
            self.levels.append(childLevel)

        for levelNodes in self.levels[::-1]:
            for node in levelNodes:
                if node.upChild != None:
                    binomialValue = (math.exp(-self.r * self.T / self.n)
                                     * (self.p * node.upChild.value + (1 - self.p) * node.downChild.value))
                else:
                    binomialValue = 0.0
                exerciseValue = self.K - node.price
                node.value = max(binomialValue, exerciseValue)

        putPrice = self.root.value
        delta = (self.root.upChild.value - self.root.downChild.value) / (self.root.upChild.price - self.root.downChild.price)

        return putPrice, delta


if __name__ == "__main__":

    ### Test Case
    # putPrice, delta = AmericanPutCRRBinomialTree(S=100., K=105., r=0.03, s=0.25, T=1, n=300).price()
    # print("put price:", putPrice)
    # print("delta:", delta)

    S = input("S:")  # spotPrice
    K = input("K:")  # strikePrice
    r = input("r:")  # riskFreeRate
    s = input("s:")  # volatility
    T = input("T:")  # years
    n = input("n:")  # numberOfPeriods

    putPrice, delta = AmericanPutCRRBinomialTree(S=float(S), K=float(K), r=float(r), s=float(s), T=float(T), n=int(n)).price()
    print("put price:", putPrice)
    print("delta:", delta)