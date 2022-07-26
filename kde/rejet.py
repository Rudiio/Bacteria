import numpy as np
import matplotlib.pyplot as plt
import joblib
from scipy import stats
import random

# Distribution à fit
incre = joblib.load("./increment.pkl")
x = np.linspace(0,10,100)
y = incre(x)
indmax = 32
mu = 3.23
# Gaussienne

mu=3.23
sigma = 1.4

L = []

for i in range(10000):
    L.append(random.gauss(mu,sigma))

c=2

def g(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)

def rejet():
    mu=3.23
    sigma=1.4
    c=2
    t = False
    while not t:
        x = random.gauss(mu,sigma)
        u = random.random()

        if(c*g(x,mu,sigma)*u < incre(x)):
            t = True
            return x 

test  = np.zeros(10000)

for i in range(10000):
    # print(i)
    test[i] = rejet()

# gauss = stats.gaussian_kde(L)
plt.hist(test,bins=100,density=True)
plt.plot(x,y,label="increment")
# plt.plot(x,c*gauss(x),label="gauss")
plt.legend()
plt.show()