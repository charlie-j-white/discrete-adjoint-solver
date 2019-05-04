import numpy as np
import matplotlib.pyplot as plt


# define column extraction 
def column(matrix, i):
        return [row[i] for row in matrix]


# read in file
solution = np.loadtxt('sensitivities.plt', skiprows=1)


# split into variables
X = column(solution, 0)
FD = column(solution, 1)
ADJ = column(solution, 2)
ERR = column(solution, 3)


# plot stuff
plt.title('calculated sensitivities')
a, = plt.plot(X,FD, '-', label='FD')
b, = plt.plot(X,ADJ, '-', label='ADJ')
#c, = plt.plot(X,ERR, '-', label='ERR')
plt.legend(handles=[a,b])
plt.show()



