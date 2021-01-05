import matplotlib.pyplot as plt

with open('displacement.txt', 'r') as f:
    for Line in f:
        l = Line.strip().split(",")
dsp = [float(l[i]) for i in range(len(l) - 1)] 
print(dsp)
with open('grid.txt', 'r') as f:
    for Line in f:
        l = Line.strip().split(",")
grid = [float(l[i]) for i in range(len(l) - 1)] 
print(dsp)

with open('Moment.txt', 'r') as f:
    for Line in f:
        l = Line.strip().split(",")
M = [float(l[i]) for i in range(len(l) - 1)] 
print(grid)

with open('Forces.txt', 'r') as f:
    for Line in f:
        l = Line.strip().split(",")
N = [float(l[i]) for i in range(len(l) - 1)] 
print(grid)



plt.figure(figsize=(10,10))
plt.plot(grid, dsp, "-*r")
plt.grid()
plt.show()

plt.figure(figsize=(10,10))
plt.plot(M, "-b")
plt.grid()
plt.show()

plt.figure(figsize=(10,10))
plt.plot(N, "-m")
plt.grid()
plt.show()