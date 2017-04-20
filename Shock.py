import numpy as np
import matplotlib.pyplot as plt
u=np.loadtxt("ushock.dat")
p=np.loadtxt("pshock.dat")
rho=np.loadtxt("rhoshock.dat")
x=np.linspace(0,1,len(u))
plt.figure(figsize=(20,20))
plt.plot(x,u)
plt.savefig("u.pdf")
plt.close()
plt.figure(figsize=(20,20))
plt.plot(x,p)
plt.savefig("p.pdf")
plt.close()
plt.figure(figsize=(20,20))
plt.plot(x,rho)
plt.savefig("rho.pdf")
plt.close()

