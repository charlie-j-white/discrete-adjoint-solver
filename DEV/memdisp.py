import numpy as np
import matplotlib.pyplot as plt


# define column extraction 
def column(matrix, i):
        return [row[i] for row in matrix]


# read in file
solution = np.loadtxt('memory.plt', skiprows=2)


# split into variables
Index = column(solution, 0)
Time = column(solution, 1)
Mem = column(solution, 2)


# plot stuff
plt.title('calculated sensitivities')
a, = plt.plot(Index,Mem, '-', label='Memory')
plt.legend(handles=[a])
plt.show()



