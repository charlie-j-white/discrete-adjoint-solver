import numpy as np
import matplotlib.pyplot as plt


# define column extraction 
def column(matrix, i):
        return [row[i] for row in matrix]


# read in file
solution = np.loadtxt('transect.plt', skiprows=1)


# split into variables
X = column(solution, 0)
u1 = column(solution, 1)
u2 = column(solution, 2)
pres = column(solution, 3)


# plot stuff
plt.title('compressible flow in a duct')
a, = plt.plot(X,u1, '-',label='u1')
b, = plt.plot(X,u2, '-',label='u2')
c, = plt.plot(X,pres, '-',label='pres')
plt.legend(handles=[a,b,c])
plt.show()



