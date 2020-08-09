import math
from scipy import special
import matplotlib.pyplot as plt
import numpy as np

N = 10
L = 1
r = 0.001
def K(s):
    return np.exp(-np.power(s, 2)/2)/np.power(s, 2)

def D(n, s):
    if n == 0:
        return math.sqrt(2*math.pi)*(1-special.erf(s/math.sqrt(2)))/(2*s)
    if n == 1:
        return K(s)
    else:
        return (n-1)*D(n-2, s)/np.power(s, 2) + K(s)

s = np.arange(0.1, L, r)
f = D(N, s)
#plt.plot(s, f)
#plt.show()
