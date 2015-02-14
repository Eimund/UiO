import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pow, pi
from scipy.stats import norm, rayleigh

def FromVMDFile(filename) :
    fil = open(filename, 'r')
    line = fil.read().split('\n')
    x = [[0 for i in range(3)] for i in range(len(line)-2)]
    v = [[0 for i in range(3)] for i in range(len(line)-2)]
    for i in range(len(line)-2) :
        buf = line[i+2].split(' ')
        for j in range(3) :
            x[i][j] = float(buf[j+1])
            v[i][j] = float(buf[j+4])     
    fil.close()
    return x,v

x,v = FromVMDFile('build-project1-Desktop-Debug/b.xyz')
mu = [0 for i in range(4)]
sigma = [0 for i in range(4)]
vel = [0 for i in range(len(v))]  
for i in range(len(x)) :
    for j in range(3) :
        mu[j] += v[i][j]
        vel[i] += v[i][j]*v[i][j]
    vel[i] = sqrt(vel[i])
    mu[3] += vel[i]
for j in range(4) :
    mu[j] /= len(v)
for i in range(len(x)) :
    for j in range(3) :
        sigma[j] += pow(v[i][j]-mu[j],2)
    sigma[3] += pow(vel[i]-mu[3],2)
for j in range(3) :
    sigma[j] = sqrt(sigma[j]/len(v))
sigma[3] = sqrt(2/pi)*mu[3]

range = np.arange(-5, 5, 0.01)
plt.figure(1)
plt.plot(range, norm.pdf(range, mu[0], sigma[0]))
plt.xlabel('$v_x$ [angstrom/ps]')
plt.ylabel('probability')
plt.title('$\mu$=%g, $\sigma$=%g'%(mu[0],sigma[0]))
plt.figure(2)
plt.plot(range, norm.pdf(range, mu[1], sigma[1]))
plt.xlabel('$v_y$ [angstrom/ps]')
plt.ylabel('probability')
plt.title('$\mu$=%g, $\sigma$=%g'%(mu[1],sigma[1]))
plt.figure(3)
plt.plot(range, norm.pdf(range, mu[2], sigma[2]))
plt.xlabel('$v_z$ [angstrom/ps]')
plt.ylabel('probability')
plt.title('$\mu$=%g, $\sigma$=%g'%(mu[2],sigma[2]))
plt.figure(4)
range = np.arange(0, 6, 0.01)
plt.plot(range, rayleigh.pdf(range, scale=sigma[3]))
plt.xlabel('$v$ [angstrom/ps]')
plt.ylabel('probability')
plt.title('mean=%g, $\sigma$=%g'%(mu[3],sigma[3]))
plt.show()
