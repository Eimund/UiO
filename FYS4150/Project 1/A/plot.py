import matplotlib.pyplot as plt
f = open('build-project1-Desktop-Debug/result_1000.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(1002)] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
t = data[0]
d = data[1]
plt.figure(1)
plt.plot(t,d)
for i in range(2,len(line)) :
    plt.plot(data[0],data[i])
plt.title('n=1000')
plt.xlabel('x')
plt.ylabel('u')
plt.legend(['exact','4n'])
f.close()
f = open('build-project1-Desktop-Debug/result_100.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(102)] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
plt.figure(2)
plt.plot(t,d)
for i in range(2,len(line)) :
    plt.plot(data[0],data[i])
plt.title('n=100')
plt.xlabel('x')
plt.ylabel('u')
plt.legend(['exact','4n'])
f.close()
f = open('build-project1-Desktop-Debug/result_10.dat', 'r')
line = f.read().split('\n')
data=[[0 for i in range(12)] for i in range(len(line))]
for i in range(len(line)) :
    buf = line[i].split('\t')
    for j in range(len(buf)) :
        data[i][j] = float(buf[j])
plt.figure(3)
plt.plot(t,d)
for i in range(2,len(line)) :
    plt.plot(data[0],data[i])
plt.title('n=10')
plt.xlabel('x')
plt.ylabel('u')
plt.legend(['exact','4n'])
f.close()
plt.show()
