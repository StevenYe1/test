import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
from Dn_plots import D

N = 2
N1 = N-1
N2 = N-1
N3 = N-1
P = 1
const = 1 #the value of e*hbar/2m
#coefficients
index = (N+1, N, 2*N+1, N1+1, N2+1, N3+1, N+1, N)
c = np.zeros(index, dtype = np.complex)
#exponents
index2 = (N+1, N)
a = np.zeros(index2)
#normalization factors
norm = np.zeros(index2)
#first three indices denote alpha (occupation number), the next three denote the spacial moments, the last two denotes the exponent
#initialize expansion coefficients given on page 235
a[1][0] = 8.80
a[2][0] = 0.252
a[2][1] = 0.385
for n in range(1, N+1):
    for l in range(n):
        norm[n, l] = (2*np.power(2*a[n][l], 3/4)/np.power(np.pi, 1/4))*np.sqrt(np.power(2, l)/sp.factorial2(2*l+1))*np.power(np.sqrt(2*a[n][l]), l)
c[2][1][2][1][0][0][2][1] = norm[2][1]*np.sqrt(3/8*np.pi)
c[2][1][2][0][1][0][2][1] = norm[2][1]*np.sqrt(3/8*np.pi)*1j
def m(i, x, y, z):
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
                                                        if (a[p1][q1] + a[p2][q2]) != 0:
                                                            if i == 1:
                                                                result = result - y*2*np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*(a[p1][q1]-a[p2][q2])*np.power(x, n1+m1)*np.power(y, n2+m2)*np.power(z, n3+m3+1)*D(n1+n2+n3+m1+m2+m3+2, np.sqrt(2*(a[n][l]+a[n][l])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
                                                                result = result + z*2*np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*(a[p1][q1]-a[p2][q2])*np.power(x, n1+m1)*np.power(y, n2+m2+1)*np.power(z, n3+m3)*D(n1+n2+n3+m1+m2+m3+2, np.sqrt(2*(a[n][l]+a[n][l])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
                                                                if (n3 + m3) > 0:
                                                                    result = result + y*np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*(m3 - n3)*np.power(x, n1+m1)*np.power(y, n2+m2)*np.power(z, n3+m3-1)*D(n1+n2+n3+m1+m2+m3, np.sqrt(2*(a[p1][q1]+a[p2][q2])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
                                                                if (n2 + m2) > 0:
                                                                    result = result - z*np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*(m2 - n2)*np.power(x, n1+m1)*np.power(y, n2+m2-1)*np.power(z, n3+m3)*D(n1+n2+n3+m1+m2+m3, np.sqrt(2*(a[p1][q1]+a[p2][q2])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
                                                            if i == 2:
                                                                result = result - z*2*np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*(a[p1][q1]-a[p2][q2])*np.power(x, n1+m1+1)*np.power(y, n2+m2)*np.power(z, n3+m3)*D(n1+n2+n3+m1+m2+m3+2, np.sqrt(2*(a[n][l]+a[n][l])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
                                                                result = result + x*2*np.conj(c[n][l][m][n1][n2][n3][p1][q1])*c[n][l][m][m1][m2][m3][p2][q2]*(a[p1][q1]-a[p2][q2])*np.power(x, n1+m1)*np.power(y, n2+m2)*np.power(z, n3+m3+1)*D(n1+n2+n3+m1+m2+m3+2, np.sqrt(2*(a[n][l]+a[n][l])*(np.power(x, 2)+np.power(y, 2)+np.power(z, 2))))
                                                                if (n1 + m1) > 0:
                                                                    result = result + z * np.conj(c[n][l][m][n1][n2][n3][p1][q1]) * c[n][l][m][m1][m2][m3][p2][q2] * (m1 - n1) * np.power(x, n1 + m1-1) * np.power(y, n2 + m2) * np.power(z, n3 + m3) * D(n1 + n2 + n3 + m1 + m2 + m3, np.sqrt(2 * (a[p1][q1] + a[p2][q2]) * (np.power(x, 2) + np.power(y, 2) + np.power(z, 2))))
                                                                if (n3 + m3) > 0:
                                                                    result = result - x * np.conj(c[n][l][m][n1][n2][n3][p1][q1]) * c[n][l][m][m1][m2][m3][p2][q2] * (m3 - n3) * np.power(x, n1 + m1) * np.power(y, n2 + m2) * np.power(z, n3 + m3-1) * D(n1 + n2 + n3 + m1 + m2 + m3, np.sqrt(2 * (a[p1][q1] + a[p2][q2]) * (np.power(x, 2) + np.power(y, 2) + np.power(z, 2))))
                                                            if i == 3:
                                                                result = result - x * 2 * np.conj(c[n][l][m][n1][n2][n3][p1][q1]) * c[n][l][m][m1][m2][m3][p2][q2] * (a[p1][q1] - a[p2][q2]) * np.power(x, n1 + m1) * np.power(y, n2 + m2+1) * np.power(z, n3 + m3) * D(n1 + n2 + n3 + m1 + m2 + m3 + 2, np.sqrt(2 * (a[n][l] + a[n][l]) * (np.power(x, 2) + np.power(y, 2) + np.power(z, 2))))
                                                                result = result + y * 2 * np.conj(c[n][l][m][n1][n2][n3][p1][q1]) * c[n][l][m][m1][m2][m3][p2][q2] * (a[p1][q1] - a[p2][q2]) * np.power(x, n1 + m1+1) * np.power(y, n2 + m2) * np.power(z, n3 + m3) * D(n1 + n2 + n3 + m1 + m2 + m3 + 2, np.sqrt(2 * (a[n][l] + a[n][l]) * (np.power(x, 2) + np.power(y, 2) + np.power(z, 2))))
                                                                if (n2 + m2) > 0:
                                                                    result = result + x * np.conj(c[n][l][m][n1][n2][n3][p1][q1]) * c[n][l][m][m1][m2][m3][p2][q2] * (m2 - n2) * np.power(x, n1 + m1) * np.power(y, n2 + m2-1) * np.power(z, n3 + m3) * D(n1 + n2 + n3 + m1 + m2 + m3, np.sqrt(2 * (a[p1][q1] + a[p2][q2]) * (np.power(x, 2) + np.power(y, 2) + np.power(z, 2))))
                                                                if (n1 + m1) > 0:
                                                                    result = result - y * np.conj(c[n][l][m][n1][n2][n3][p1][q1]) * c[n][l][m][m1][m2][m3][p2][q2] * (m1 - n1) * np.power(x, n1 + m1-1) * np.power(y, n2 + m2) * np.power(z, n3 + m3) * D(n1 + n2 + n3 + m1 + m2 + m3, np.sqrt(2 * (a[p1][q1] + a[p2][q2]) * (np.power(x, 2) + np.power(y, 2) + np.power(z, 2))))
    result = result * (-1j) * const
    return np.real(result)


fig = plt.figure(1)
ax = fig.gca(projection = '3d')
z = np.zeros(11*81-1)
x = np.zeros(11*81-1)
y = np.zeros(11*81-1)
ind = 0
scale = 0.4
for k in range(11):
    z[ind] = scale*(k-5)
    x[ind] = 0
    y[ind] = 0
    ind = ind + 1
    if k == 5:
        ind = ind - 1
    for i in range(5):
        for j in range(8*i):
            x1 = scale*i*np.cos(2*j*np.pi/(8*i))
            y1 = scale*i*np.sin(2*j*np.pi/(8*i))
            z[ind] = scale*(k-5)
            x[ind] = x1
            y[ind] = y1
            ind = ind + 1

u = m(1, x, y, z)
v = m(2, x, y, z)
w = m(3, x, y, z)
t = np.power(np.power(np.abs(u), 2)+np.power(np.abs(v), 2)+np.power(np.abs(w), 2), 0.25)
u = u/t
v = v/t
w = w/t
ax.quiver(x, y, z, u, v, w, length = 0.1, normalize = False)
plt.show()