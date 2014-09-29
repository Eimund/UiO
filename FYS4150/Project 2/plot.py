import numpy as np
import matplotlib.pyplot as plt
from numpy import linspace
from math import sqrt, pi, exp
f = open('build-project2-Desktop-Debug/QRw0_2_4.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
f.close()
plt.figure(1)
plt.plot(data[0],data[1])
plt.figure(2)
plt.plot(data[0],data[2])
f = open('build-project2-Desktop-Debug/QRw1_1_0.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
f.close()
plt.figure(1)
plt.plot(data[0],data[1])
plt.figure(2)
plt.plot(data[0],data[2])
f = open('build-project2-Desktop-Debug/QRw2_1_0.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
f.close()
plt.figure(1)
plt.plot(data[0],data[1])
plt.figure(2)
plt.plot(data[0],data[2])
f = open('build-project2-Desktop-Debug/QRw3_1_0.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
f.close()
plt.figure(1)
plt.plot(data[0],data[1])
plt.figure(2)
plt.plot(data[0],data[2])
plt.figure(1)
plt.xlim([0,30])
plt.title('$tqli (lib.cpp)$')
plt.xlabel('$\\rho$')
plt.ylabel('|R|')
plt.legend(['$\omega_r = 0.01$', '$\omega_r = 0.5$' ,'$\omega_r = 1$' , '$\omega_r = 5$'])
plt.figure(2)
plt.xlim([0,30])
plt.title('$QR$')
plt.xlabel('$\\rho$')
plt.ylabel('|R|')
plt.legend(['$\omega_r = 0.01$', '$\omega_r = 0.5$' ,'$\omega_r = 1$' , '$\omega_r = 5$'])

plt.figure(3)
omega = 0.01;
r = linspace(0,30,1000)
d = (sqrt(3)*omega/pi)**(0.25)*np.exp(-0.5*sqrt(3)*omega*(np.array(r)-(2*omega**2)**(-1.0/3.0))**2)
d = d/sum(d*d)
plt.plot(r,d)
omega = 0.5;
d = (sqrt(3)*omega/pi)**(0.25)*np.exp(-0.5*sqrt(3)*omega*(np.array(r)-(2*omega**2)**(-1.0/3.0))**2)
d = d/sum(d*d)
plt.plot(r,d)
omega = 1;
d = (sqrt(3)*omega/pi)**(0.25)*np.exp(-0.5*sqrt(3)*omega*(np.array(r)-(2*omega**2)**(-1.0/3.0))**2)
d = d/sum(d*d)
plt.plot(r,d)
omega = 5;
d = (sqrt(3)*omega/pi)**(0.25)*np.exp(-0.5*sqrt(3)*omega*(np.array(r)-(2*omega**2)**(-1.0/3.0))**2)
d = d/sum(d*d)
plt.plot(r,d)
plt.title('$Exact$')
plt.xlabel('$\\rho$')
plt.ylabel('|R|')
plt.legend(['$\omega_r = 0.01$', '$\omega_r = 0.5$' ,'$\omega_r = 1$' , '$\omega_r = 5$'])
plt.show()