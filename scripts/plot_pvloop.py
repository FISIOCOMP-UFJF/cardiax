import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
data = np.loadtxt(filename)

plt.plot(data[:,2], data[:,1])

# plt.figure(0)
# plt.plot(data[:,0], data[:,1], label='Plv')
# plt.plot(data[:,0], data[:,4], label='Part')
# plt.plot(data[:,0], data[:,5], label='Pven')
# plt.plot(data[:,0], data[:,6], label='PLA')

# plt.legend()

# plt.figure(1)
# plt.plot(data[:,0], data[:,2], label='Vlv')

# plt.figure(2)
# plt.plot(data[:,0], data[:,10], label='qao')
# plt.plot(data[:,0], data[:,11], label='qmv')
# plt.plot(data[:,0], data[:,12], label='qven')
# plt.plot(data[:,0], data[:,13], label='qper')
# plt.legend()
# plt.figure(3)
# plt.plot(data[:,0], data[:,3], label='Ta')
# plt.legend()
# plt.figure(4)
# plt.plot(data[:,0], data[:,-1], label='C')
# plt.legend()
#np.savetxt('elastancia.txt', np.array([data[:,0], data[:,-1]]).T)
plt.show()
