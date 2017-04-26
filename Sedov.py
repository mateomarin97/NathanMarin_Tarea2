import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rho=np.loadtxt("rhosedov.dat")
x=np.linspace(0,128,len(rho))
plt.figure(figsize=(20,20))
plt.plot(x,rho)
plt.savefig("rhosedov.pdf")
plt.closefig()
