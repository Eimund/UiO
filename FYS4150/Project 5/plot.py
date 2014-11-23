import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as nu

f = open('build-project5-Desktop-Debug/Exact_2D.dat', 'r')
line = f.read().split('\n')
buf = line[0].split('\t')

x1=[0 for i in range(len(line[0].split('\t')))]
for i in range(len(buf)) :
    x1[i] = float(buf[i])
buf = line[1].split('\t')
x2=[0 for i in range(len(line[1].split('\t')))]
for i in range(len(buf)) :
    x2[i] = float(buf[i])
data=[[0 for i in range(len(line[2].split('\t')))] for i in range(len(line)-2)] 
for i in range(len(line)-2) :
    buf = line[i+2].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
        
f.close()
        
x, y = nu.meshgrid(x2,x1) 

fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(x,y, data, rstride=4, cstride=4)

f = open('build-project5-Desktop-Debug/Explicit_2D.dat', 'r')
line = f.read().split('\n')
buf = line[0].split('\t')

x1=[0 for i in range(len(line[0].split('\t')))]
for i in range(len(buf)) :
    x1[i] = float(buf[i])
buf = line[1].split('\t')
x2=[0 for i in range(len(line[1].split('\t')))]
for i in range(len(buf)) :
    x2[i] = float(buf[i])
data=[[0 for i in range(len(line[2].split('\t')))] for i in range(len(line)-2)] 
for i in range(len(line)-2) :
    buf = line[i+2].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
        
f.close()
        
x, y = nu.meshgrid(x2,x1) 

fig = plt.figure(2)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(x,y, data, rstride=4, cstride=4)

f = open('build-project5-Desktop-Debug/Jacobi_2D.dat', 'r')
line = f.read().split('\n')
buf = line[0].split('\t')

x1=[0 for i in range(len(line[0].split('\t')))]
for i in range(len(buf)) :
    x1[i] = float(buf[i])
buf = line[1].split('\t')
x2=[0 for i in range(len(line[1].split('\t')))]
for i in range(len(buf)) :
    x2[i] = float(buf[i])
data=[[0 for i in range(len(line[2].split('\t')))] for i in range(len(line)-2)] 
for i in range(len(line)-2) :
    buf = line[i+2].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
        
f.close()
        
x, y = nu.meshgrid(x2,x1) 

fig = plt.figure(3)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(x,y, data, rstride=4, cstride=4)

plt.show()
         
