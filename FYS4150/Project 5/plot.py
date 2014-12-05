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

x, data = FromFile('build-project5-Desktop-Debug/Exact_2D.dat')
X, Y = nu.meshgrid(x[1],x[0]) 
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(X,Y, data, rstride=4, cstride=4)

x, data = FromFile('build-project5-Desktop-Debug/Explicit_2D.dat')
X, Y = nu.meshgrid(x[1],x[0]) 
fig = plt.figure(2)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(X,Y, data, rstride=4, cstride=4)

x, data = FromFile('build-project5-Desktop-Debug/Jacobi_2D.dat')
X, Y = nu.meshgrid(x[1],x[0]) 
fig = plt.figure(3)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(X,Y, data, rstride=4, cstride=4)

x, data = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Uniform.dat')
fig = plt.figure(4)
plt.plot(x[0],data)

x, data = FromFile('build-project5-Desktop-Debug/MonteCarlo_1D_Gaussian.dat')
fig = plt.figure(5)
plt.plot(x[0],data)

x, data = FromFile('build-project5-Desktop-Debug/MonteCarlo_2D_Gaussian.dat')
X, Y = nu.meshgrid(x[1],x[0]) 
fig = plt.figure(6)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(X,Y, data, rstride=4, cstride=4)

x, data = FromFile('build-project5-Desktop-Debug/Metropolis_1D.dat')
fig = plt.figure(7)
plt.plot(x[0],data)

x, data = FromFile('build-project5-Desktop-Debug/Exact_1D.dat')
fig = plt.figure(8)
plt.plot(x[0],data)

x, data = FromFile('build-project5-Desktop-Debug/Metropolis_2D.dat')
X, Y = nu.meshgrid(x[1],x[0]) 
fig = plt.figure(9)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(X,Y, data, rstride=4, cstride=4)

plt.show()
         
