import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rho3=np.loadtxt("rho3sedov.dat")
#P3=np.loadtxt("p3sedov.dat")
#u3=np.loadtxt("u3sedov.dat")
rho2=np.loadtxt("rho2sedov.dat")
#P2=np.loadtxt("p2sedov.dat")
#u2=np.loadtxt("u2sedov.dat")
rho1=np.loadtxt("rho1sedov.dat")
#P1=np.loadtxt("p1sedov.dat")
#u1=np.loadtxt("u1sedov.dat")
x=np.linspace(0,128,len(rho3))
plt.figure(figsize=(20,20))
plt.plot(x,rho3)
plt.xlim(0,128)
plt.savefig("rho120sedov.pdf")
plt.close()

plt.figure(figsize=(20,20))
plt.plot(x,rho2)
plt.xlim(0,128)
plt.savefig("rho60sedov.pdf")
plt.close()

plt.figure(figsize=(20,20))
plt.plot(x,rho1)
plt.xlim(0,128)
plt.savefig("rho10sedov.pdf")
plt.close()


