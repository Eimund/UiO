import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as nu

def Space(n) :
    space = [0 for i in range(n[0])]
    if(len(n) > 1) :
        for i in range(n[0]) :
            space[i] = Space(n[1:])
    return space

def ReadLine(line, i, n) :
    space = Space([n[0]])
    if(len(n) > 1) :
        for j in range(n[0]) :
            space[j], i = ReadLine(line, i, n[1:])
    else :
        buf = line[i].split('\t')
        for j in range(n[0]) :
            space[j] = float(buf[j])
        i = i+1
    return space, i    

def FromFile(filename) :
    fil = open(filename, 'r')
    line = fil.read().split('\n')
    n = Space([int(line[0])])
    x = Space([len(n)])
    for i in range(len(n)) :
        buf = line[i+1].split('\t')
        n[i] = len(buf)
        x[i] = Space([n[i]])
        for j in range(n[i]) :
            x[i][j] = float(buf[j])
    space, i = ReadLine(line[len(n)+1:], 0, n)      
    fil.close()
    return x, space

x1, data1 = FromFile('build-project5-Desktop-Debug/Exact_1D_t0.01.dat')
x2, data2 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform_t0.01_dt1e-3_n1.dat')
x3, data3 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform_t0.01_dt1e-3_n100.dat')
x4, data4 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform_t0.01_dt1e-3_n10000.dat')
fig = plt.figure(1)
plt.subplot(1,3,1)
plt.plot(x1[0],data1,'c',x2[0],data2,'g', x3[0], data3,'b', x4[0], data4,'r')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Monte Carlo Uniform: t=0.01, dt=1e-3')
plt.ylim([0,1])
plt.legend(['Exact', '1 experiment', '100 experiments', '10000 experiments'])
x2, data2 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform_t0.01_dt1e-5_n1.dat')
x3, data3 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform_t0.01_dt1e-5_n100.dat')
x4, data4 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform_t0.01_dt1e-5_n10000.dat')
plt.subplot(1,3,2)
plt.plot(x1[0],data1,'c',x2[0],data2,'g', x3[0], data3,'b', x4[0], data4,'r')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Monte Carlo Uniform: t=0.01, dt=1e-5')
plt.ylim([0,1])
plt.legend(['Exact', '1 experiment', '100 experiments', '10000 experiments'])
x2, data2 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian_t0.01_dt1e-4_n1.dat')
x3, data3 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian_t0.01_dt1e-4_n100.dat')
x4, data4 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian_t0.01_dt1e-4_n10000.dat')
plt.subplot(1,3,3)
plt.plot(x1[0],data1,'c',x2[0],data2,'g', x3[0], data3,'b', x4[0], data4,'r')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Monte Carlo Gaussian: t=0.01, dt=1e-4')
plt.ylim([0,1])
plt.legend(['Exact', '1 experiment', '100 experiments', '10000 experiments'])

x1, data1 = FromFile('build-project5-Desktop-Debug/Exact_1D_t1.dat')
x2, data2 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian_t1_dt1e-3_n1.dat')
x3, data3 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian_t100_dt1e-3_n1.dat')
x4, data4 = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian_t1000_dt1e-3_n1.dat')
fig = plt.figure(2)
plt.subplot(1,3,1)
plt.plot(x1[0],data1,'c',x2[0],data2,'g',x3[0], data3,'b', x4[0], data4,'r')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Monte Carlo Gaussian: t=1, dt=1e-3')
plt.ylim([0,1])
plt.legend(['Exact', '1 experiment', '100 experiments','1000 experiments'])

#x, data = FromFile('build-project5-Desktop-Debug/Exact_2D.dat')
#X, Y = nu.meshgrid(x[1],x[0]) 
#fig = plt.figure(1)
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.plot_surface(X,Y, data, rstride=4, cstride=4)

#x, data = FromFile('build-project5-Desktop-Debug/Explicit_2D.dat')
#X, Y = nu.meshgrid(x[1],x[0]) 
#fig = plt.figure(2)
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.plot_surface(X,Y, data, rstride=4, cstride=4)

#x, data = FromFile('build-project5-Desktop-Debug/Jacobi_2D.dat')
#X, Y = nu.meshgrid(x[1],x[0]) 
#fig = plt.figure(3)
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.plot_surface(X,Y, data, rstride=4, cstride=4)

#x, data = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform.dat')
#fig = plt.figure(4)
#plt.plot(x[0],data)

#x, data = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian.dat')
#fig = plt.figure(5)
#plt.plot(x[0],data)

#x, data = FromFile('build-project5-Desktop-Debug/MonteCarlo_2D_Gaussian.dat')
#X, Y = nu.meshgrid(x[1],x[0]) 
#fig = plt.figure(6)
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.plot_surface(X,Y, data, rstride=4, cstride=4)

#x, data = FromFile('build-project5-Desktop-Debug/Metropolis_1D.dat')
#fig = plt.figure(7)
#plt.plot(x[0],data)

#x, data = FromFile('build-project5-Desktop-Debug/Exact_1D.dat')
#fig = plt.figure(8)
#plt.plot(x[0],data)

#x, data = FromFile('build-project5-Desktop-Debug/Metropolis_2D.dat')
#X, Y = nu.meshgrid(x[1],x[0]) 
#fig = plt.figure(9)
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.plot_surface(X,Y, data, rstride=4, cstride=4)

plt.show()
         
