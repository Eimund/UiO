import matplotlib.pyplot as plt
from _dbus_bindings import String

f = open('build-project3-Desktop-Debug/Verlet_Sun_and_Earth_10.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=1, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=10 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(num=2, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=10')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/Verlet_Sun_and_Earth_20.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(1)
plt.subplot(2,2,2)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=20 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(2)
plt.subplot(2,2,2)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=20')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/Verlet_Sun_and_Earth_50.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(1)
plt.subplot(2,2,3)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=50 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(2)
plt.subplot(2,2,3)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=50')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/Verlet_Sun_and_Earth_100.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(1)
plt.subplot(2,2,4)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=100 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(2)
plt.subplot(2,2,4)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=100')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_10.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=3, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=10 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(num=4, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=10')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_20.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=3)
plt.subplot(2,2,2)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=20 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(num=4)
plt.subplot(2,2,2)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=20')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_50.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=3)
plt.subplot(2,2,3)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=50 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(num=4)
plt.subplot(2,2,3)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=50')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_100.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=3)
plt.subplot(2,2,4)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=100 - 1yr')
plt.legend([names[0], names[1]])
plt.figure(num=4)
plt.subplot(2,2,4)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=100')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_20.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=5, figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(data[0],data[1])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Sun - RK4 n=20 - 1yr')
plt.figure(num=6, figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=20')
plt.legend(["%s - x" % names[0], "%s - y" % names[0]])

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_20.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=5, figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(data[0],data[1],'b')
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Sun and Earth - RK4 n=20 - 1yr')
plt.legend([names[0]])
plt.figure(num=6, figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Sun and Earth - RK4 n=20')
plt.legend(["%s - x" % names[0], "%s - y" % names[0]])

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_50.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=5, figsize=(10,10))
plt.subplot(2,1,2)
plt.plot(data[0],data[1])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Sun and Earth - RK4 n=50 - 1yr')
plt.legend([names[0]])
plt.figure(num=6, figsize=(10,10))
plt.subplot(2,1,2)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Sun and Earth - RK4 n=50')
plt.legend(["%s - x" % names[0], "%s - y" % names[0]])

f = open('build-project3-Desktop-Debug/RK4_Energy_momentum_100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
f.close()
plt.figure(num=7)
plt.plot(data[0],data[1], data[0],data[2], data[0],data[3], data[0],data[4])
plt.xlabel('[yr]')
plt.title('Sun and Earth - RK4 n=100')
plt.legend(['Kinetic energy','Potential energy', 'Momentum x', 'Momentum y'],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_and_Earth_escape_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=8)
plt.plot(data[4],data[6])
plt.xlabel('[AU]')
plt.ylabel('[AU/yr]')
plt.title('Escape velocity of Earth from the Sun - RK4 n=1000')

f = open('build-project3-Desktop-Debug/Verlet_Sun_Earth_Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=9, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-6,6])
plt.ylim([-6,6])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=1000 - 11.8618yr - Jupiter x1')
plt.legend([names[0], names[1], names[2]])
plt.figure(num=10, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=1000 - Jupiter x1')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=4)

f = open('build-project3-Desktop-Debug/Verlet_Sun_Earth_10Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=9)
plt.subplot(2,2,2)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-6,6])
plt.ylim([-6,6])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=1000 - 11.8618yr - Jupiter x10')
plt.legend([names[0], names[1], names[2]])
plt.figure(num=10)
plt.subplot(2,2,2)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=1000 - Jupiter x10')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=4)

f = open('build-project3-Desktop-Debug/Verlet_Sun_Earth_100Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=9)
plt.subplot(2,2,3)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-6,6])
plt.ylim([-6,6])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=1000 - 11.8618yr - Jupiter x100')
plt.legend([names[0], names[1], names[2]],loc=4)
plt.figure(num=10)
plt.subplot(2,2,3)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=1000 - Jupiter x100')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=4)

f = open('build-project3-Desktop-Debug/Verlet_Sun_Earth_1000Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=9)
plt.subplot(2,2,4)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-20,20])
plt.ylim([-20,20])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=1000 - 11.8618yr - Jupiter x1000')
plt.legend([names[0], names[1], names[2]],loc=1)
plt.figure(num=10)
plt.subplot(2,2,4)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=1000 - Jupiter x1000')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=2)

f = open('build-project3-Desktop-Debug/RK4_Sun_Earth_Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=11, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-6,6])
plt.ylim([-6,6])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=1000 - 11.8618yr - Jupiter x1')
plt.legend([names[0], names[1], names[2]])
plt.figure(num=12, figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=1000 - Jupiter x1')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_Earth_10Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=11)
plt.subplot(2,2,2)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-6,6])
plt.ylim([-6,6])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=1000 - 11.8618yr - Jupiter x10')
plt.legend([names[0], names[1], names[2]])
plt.figure(num=12)
plt.subplot(2,2,2)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=1000 - Jupiter x10')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_Earth_100Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=11)
plt.subplot(2,2,3)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-6,6])
plt.ylim([-6,6])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=1000 - 11.8618yr - Jupiter x100')
plt.legend([names[0], names[1], names[2]],loc=4)
plt.figure(num=12)
plt.subplot(2,2,3)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=1000 - Jupiter x100')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_Sun_Earth_1000Jupiter_1000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t')))]
buf = line[0].split('\t')
for i in range(len(buf)) :
    t[i] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            data[i-j-1][k] = float(buf[k])
f.close()
plt.figure(num=11)
plt.subplot(2,2,4)
plt.plot(data[0],data[1],'o')
plt.plot(data[4],data[5])
plt.plot(data[8],data[9])
plt.xlim([-20,20])
plt.ylim([-20,20])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=1000 - 11.8618yr - Jupiter x1000')
plt.legend([names[0], names[1], names[2]],loc=1)
plt.figure(num=12)
plt.subplot(2,2,4)
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[6],'--g')
plt.plot(t,data[7],'.g')
plt.plot(t,data[10],'--r')
plt.plot(t,data[11],'.r')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=1000 - Jupiter x1000')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[2], "%s - y" % names[2]],loc=2)

f = open('build-project3-Desktop-Debug/Verlet_SolarSystem_10000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t'))/10)] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t'))/10)]
buf = line[0].split('\t')
for i in range(len(buf)) :
    if((i%10) == 0):
        t[i/10] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            if((k%10) == 0):
                data[i-j-1][k/10] = float(buf[k])
f.close()
plt.figure(num=13, figsize=(20,20))
plt.plot(data[0],data[1],'o')
plt.plot(data[12],data[13])
plt.plot(data[16],data[17])
plt.plot(data[4],data[5])
plt.plot(data[20],data[21])
plt.plot(data[8],data[9])
plt.plot(data[24],data[25])
plt.plot(data[28],data[29], color = '#ff7700')
plt.plot(data[32],data[33], color = '#ff0077')
plt.plot(data[36],data[37], color = '#77ff00')
plt.xlim([-40,40])
plt.ylim([-40,40])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('Verlet n=40000 - 250yr')
plt.legend([names[0], names[3], names[4], names[1], names[5], names[2], names[6], names[7], names[8], names[9]],loc=4)
plt.figure(num=14, figsize=(20,20))
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[14],'--g')
plt.plot(t,data[15],'.g')
plt.plot(t,data[18],'--r')
plt.plot(t,data[19],'.r')
plt.plot(t,data[6],'--c')
plt.plot(t,data[7],'.c')
plt.plot(t,data[22],'--m')
plt.plot(t,data[23],'.m')
plt.plot(t,data[10],'--y')
plt.plot(t,data[11],'.y')
plt.plot(t,data[26],'--k')
plt.plot(t,data[27],'.k')
plt.plot(t,data[30], color = '#ff7700', linestyle = '--')
plt.plot(t,data[31], color = '#ff7700', marker = '.', linestyle='')
plt.plot(t,data[34], color = '#ff0077', linestyle = '--')
plt.plot(t,data[35], color = '#ff0077', marker = '.', linestyle='')
plt.plot(t,data[38], color = '#77ff00', linestyle = '--')
plt.plot(t,data[39], color = '#77ff00', marker = '.', linestyle='')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('Verlet n=40000')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[3], "%s - y" % names[3], "%s - x" % names[4], "%s - y" % names[4], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[5], "%s - y" % names[5], "%s - x" % names[2], "%s - y" % names[2], "%s - x" % names[6], "%s - y" % names[6], "%s - x" % names[7], "%s - y" % names[7], "%s - x" % names[8], "%s - y" % names[8], "%s - x" % names[9], "%s - y" % names[9]],loc=4)

f = open('build-project3-Desktop-Debug/RK4_SolarSystem_10000.dat', 'r')
line = f.read().split('\n')
names = [0 for i in range((len(line)-1)/5)]
data=[[0 for i in range(len(line[0].split('\t'))/10)] for i in range(len(names)*4)]
t = [0 for i in range(len(line[0].split('\t'))/10)]
buf = line[0].split('\t')
for i in range(len(buf)) :
    if((i%10) == 0):
        t[i/10] = float(buf[i])
for i in range(len(line)-1) :
    buf = line[i+1].split('\t')
    reminder = i % 5
    j = (i-reminder)/5
    if reminder == 0 :
        names[j] = str(buf[0])
    else :
        for k in range(len(buf)) :
            if((k%10) == 0):
                data[i-j-1][k/10] = float(buf[k])
f.close()
plt.figure(num=15, figsize=(20,20))
plt.plot(data[0],data[1],'o')
plt.plot(data[12],data[13])
plt.plot(data[16],data[17])
plt.plot(data[4],data[5])
plt.plot(data[20],data[21])
plt.plot(data[8],data[9])
plt.plot(data[24],data[25])
plt.plot(data[28],data[29], color = '#ff7700')
plt.plot(data[32],data[33], color = '#ff0077')
plt.plot(data[36],data[37], color = '#77ff00')
plt.xlim([-40,40])
plt.ylim([-40,40])
plt.xlabel('[AU]')
plt.ylabel('[AU]')
plt.title('RK4 n=40000 - 250yr')
plt.legend([names[0], names[3], names[4], names[1], names[5], names[2], names[6], names[7], names[8], names[9]],loc=4)
plt.figure(num=16, figsize=(20,20))
plt.plot(t,data[2],'--b')
plt.plot(t,data[3],'.b')
plt.plot(t,data[14],'--g')
plt.plot(t,data[15],'.g')
plt.plot(t,data[18],'--r')
plt.plot(t,data[19],'.r')
plt.plot(t,data[6],'--c')
plt.plot(t,data[7],'.c')
plt.plot(t,data[22],'--m')
plt.plot(t,data[23],'.m')
plt.plot(t,data[10],'--y')
plt.plot(t,data[11],'.y')
plt.plot(t,data[26],'--k')
plt.plot(t,data[27],'.k')
plt.plot(t,data[30], color = '#ff7700', linestyle = '--')
plt.plot(t,data[31], color = '#ff7700', marker = '.', linestyle='')
plt.plot(t,data[34], color = '#ff0077', linestyle = '--')
plt.plot(t,data[35], color = '#ff0077', marker = '.', linestyle='')
plt.plot(t,data[38], color = '#77ff00', linestyle = '--')
plt.plot(t,data[39], color = '#77ff00', marker = '.', linestyle='')
plt.xlabel('[yr]')
plt.ylabel('[AU/yr]')
plt.title('RK4 n=40000')
plt.legend(["%s - x" % names[0], "%s - y" % names[0], "%s - x" % names[3], "%s - y" % names[3], "%s - x" % names[4], "%s - y" % names[4], "%s - x" % names[1], "%s - y" % names[1], "%s - x" % names[5], "%s - y" % names[5], "%s - x" % names[2], "%s - y" % names[2], "%s - x" % names[6], "%s - y" % names[6], "%s - x" % names[7], "%s - y" % names[7], "%s - x" % names[8], "%s - y" % names[8], "%s - x" % names[9], "%s - y" % names[9]],loc=4)


plt.show()