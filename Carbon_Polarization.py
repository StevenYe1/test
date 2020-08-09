import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
from Dn_plots import D

N = 2
N1 = N-1
N2 = N-1
N3 = N-1
P = 1
#coefficients
index = (N+1, N, 2*N+1, N1+1, N2+1, N3+1, N+1, N)
c = np.zeros(index, dtype = np.complex)
#exponents
index2 = (N+1, N)
a = np.zeros(index2)
#normalization factors
norm = np.zeros(index2)
#first three indices denote alpha (occupation number), the next three denote the spacial moments
#initialize expansion coefficients given on page 235
a[1][0] = 8.80
a[2][0] = 0.252
a[2][1] = 0.385
#scaling constant
const = 0.3
for n in range(1, N+1):
    for l in range(n):
        norm[n, l] = (2*np.power(2*a[n][l], 3/4)/np.power(np.pi, 1/4))*np.sqrt(np.power(2, l)/sp.factorial2(2*l+1))*np.power(np.sqrt(2*a[n][l]), l)
c[2][1][2][1][0][0][2][1] = norm[2][1]*np.sqrt(8*np.pi/3)
c[2][1][2][0][1][0][2][1] = norm[2][1]*np.sqrt(8*np.pi/3)*1j
def p(i, x, y, z):
    result = 0
    for n in range(1, N+1):
        for l in range(n):
            for m in range(2*l+1):
                for n1 in range(N1+1):
                    for n2 in range(N2+1):
                        for n3 in range(N3+1):
                            for m1 in range(N1+1):
                                for m2 in range(N2+1):
                                    for m3 in range(N3+1):
                                        for p1 in range(1, N + 1):
                                            for q1 in range(n):
                                                for p2 in range(1, N + 1):
                                                    for q2 in range(n):
                                                        v = 0
                                                        if i == 1:
                                                            v = x
                                                        if i == 2:
                                                            v = y
                                                        if i == 3:
                                                            v = z
                                                        if (a[p1][q1] + a[p2][q2]) != 0:
                                                            result = result + np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*v*np.power(x, n1+m1)*np.power(y, n2+m2)*np.power(z, n3+m3)*D(n1+n2+n3+m1+m2+m3+2, np.sqrt(2*(a[p1][q1]+a[p2][q2])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
    result = const*result
    return np.real(result)

L = 2
r = 0.4
fig = plt.figure(1)
ax = fig.gca(projection = '3d')
x, y, z = np.meshgrid((np.arange(-L, L, r)), (np.arange(-L, L, r)), (np.arange(-L, L, r)))
u = p(1, x, y, z)
v = p(2, x, y, z)
w = p(3, x, y, z)
#set zero point to zero to avoid blow-up
u[int(L/r), int(L/r), int(L/r)] = 0.1
v[int(L/r), int(L/r), int(L/r)] = 0.1
w[int(L/r), int(L/r), int(L/r)] = 0.1
u = u*np.power(np.abs(u), -0.5)
v = v*np.power(np.abs(v), -0.5)
w = w*np.power(np.abs(w), -0.5)
u[int(L/r), int(L/r), int(L/r)] = 0
v[int(L/r), int(L/r), int(L/r)] = 0
w[int(L/r), int(L/r), int(L/r)] = 0
ax.quiver(x, y, z, u, v, w, length = 0.1, normalize = False)
fig2, ax2 = plt.subplots()
s, t = np.meshgrid((np.arange(-L, L, r)), (np.arange(-L, L, r)))
f = p(1, s, t, 0)
g = p(2, s, t, 0)
f[int(L/r), int(L/r)] = 0.1
g[int(L/r), int(L/r)] = 0.1
f = f*np.power(np.abs(f), -0.5)
g = g*np.power(np.abs(g), -0.5)
f[int(L/r), int(L/r)] = 0
g[int(L/r), int(L/r)] = 0
ax2.quiver(s, t, f, g)
plt.show()