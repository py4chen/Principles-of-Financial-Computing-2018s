import math

class BinomialTree:
    class Node:

        def __init__(self, price, upParent=None, downParent=None, downChild=None, upChild=None):
            self.price = price
            self.value = None
            self.upParent = upParent
            self.downParent = downParent
            self.downChild = downChild
            self.upChild = upChild

        def insertDownChild(self, price):
            self.downChild = BinomialTree.Node(price=price, upParent=self)
            return self.downChild

        def insertUpChild(self, price):
            self.upChild = BinomialTree.Node(price=price, upParent=self)
            return self.upChild

        def insertDownExistChild(self, child):
            self.downChild = child
            child.upParent = self
            return child

        def insertUpExistChild(self, child):
            self.upChild = child
            child.downParent = self
            return child

    def __init__(self, spotPrice, upFactor, downFactor, numPeriods):
        self.root = BinomialTree.Node(spotPrice)
        self.levels = [[self.root]]

        for period in range(numPeriods):
            currentLevel = self.levels[-1]
            childLevel = [BinomialTree.Node(price = currentLevel[0].price * upFactor)]

            for currentNode in currentLevel:
                currentNode.insertUpExistChild(childLevel[-1])
                childLevel.append(currentNode.insertDownChild(currentNode.price * downFactor))

            self.levels.append(childLevel)

    def traverse(self, levelNodes, K, r, dt, p,):
        for node in levelNodes:
            assert(node.value == None)

            if node.upChild != None:
                binomialValue = (math.exp(-r * dt) * (p * node.upChild.value + (1 - p) * node.downChild.value))
            else:
                binomialValue = 0.0

            exerciseValue = K - node.price

            node.value = max(binomialValue, exerciseValue)

    def reversedLevelTraverse(self, K, r, dt, p):
        for levelNodes in self.levels[::-1]:
            self.traverse(levelNodes, K, r, dt, p)


def calAmericanPut(S, K, r, s, T, n):



    dt = T / n
    u = math.exp(s * math.sqrt(dt)) #Page 284
    d = 1 / u
    p = (math.exp(r * dt) - d) / (u - d)

    tree = BinomialTree(spotPrice=S, upFactor=u, downFactor=d, numPeriods=n)

    tree.reversedLevelTraverse(K, r, dt, p)

    firstChildren = [tree.root.upChild, tree.root.downChild]


    putPrice = tree.root.value
    delta = (firstChildren[0].value - firstChildren[1].value) / (firstChildren[0].price - firstChildren[1].price)


    print("American put price:", putPrice)
    print("Delta for the put:", delta)


if __name__ == "__main__":
    calAmericanPut(S=100., K=105., r=0.05, s=0.3, T=0.75, n=300)
    calAmericanPut(S=100., K=105., r=0.03, s=0.25, T=1, n=300)
