import matplotlib.pyplot as plt

f = open('build-project4-Desktop-Debug/Exact_t0.01_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/FE_t0.01_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/BE_t0.01_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/CN_t0.01_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
plt.ylim([0,1])
plt.xlabel('$x$')
plt.ylabel('$u$')
plt.title('$n_x = 10$, $\\alpha=0.49$, $t=0.01$')
plt.legend(['Exact','FE','BE','CN'],loc=1)

plt.figure(2)
f = open('build-project4-Desktop-Debug/Exact_t1_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/FE_t1_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/BE_t1_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/CN_t1_a0.49_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
plt.ylim([0,1])
plt.xlabel('$x$')
plt.ylabel('$u$')
plt.title('$n_x = 10$, $\\alpha=0.49$, $t=1$')
plt.legend(['Exact','FE','BE','CN'],loc=1)

plt.figure(3)
f = open('build-project4-Desktop-Debug/Exact_t0.01_a0.1_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/FE_t0.01_a0.1_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/BE_t0.01_a0.1_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/CN_t0.01_a0.1_N10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
plt.ylim([0,1])
plt.xlabel('$x$')
plt.ylabel('$u$')
plt.title('$n_x = 10$, $\\alpha=0.1$, $t=0.01$')
plt.legend(['Exact','FE','BE','CN'],loc=1)

plt.figure(4)
f = open('build-project4-Desktop-Debug/Exact_t0.01_a0.49_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/FE_t0.01_a0.49_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/BE_t0.01_a0.49_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/CN_t0.01_a0.49_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
plt.ylim([0,1])
plt.xlabel('$x$')
plt.ylabel('$u$')
plt.title('$n_x = 100$, $\\alpha=0.49$, $t=0.01$')
plt.legend(['Exact','FE','BE','CN'],loc=1)

plt.figure(5)
f = open('build-project4-Desktop-Debug/Exact_t0.01_a0.51_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/FE_t0.01_a0.51_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/BE_t0.01_a0.51_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
f = open('build-project4-Desktop-Debug/CN_t0.01_a0.51_N100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(len(line[0].split('\t')))] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])

f.close()
plt.plot(data[0],data[1])
plt.ylim([0,1])
plt.xlabel('$x$')
plt.ylabel('$u$')
plt.title('$n_x = 100$, $\\alpha=0.51$, $t=0.01$')
plt.legend(['Exact','FE','BE','CN'],loc=1)

plt.show()