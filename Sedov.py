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
dt3=np.loadtxt("t3sedov.dat")
dt2=np.loadtxt("t2sedov.dat")
dt1=np.loadtxt("t1sedov.dat")
dt3=dt3/np.sqrt(1020200)
dt2=dt2/np.sqrt(1020200)
dt1=dt1/np.sqrt(1020200)

plt.figure(figsize=(20,20))
plt.plot(x,rho3)
plt.text(20.0,1.1,"t="+str(dt3)+"s",va='center',ha='left',size=20)
plt.xlim(0,128)
plt.xlabel(r"$\mathrm{x}$ $(m)$",size=20)
plt.ylabel(r"$\mathrm{\rho}$ $\mathrm{\frac{Kg}{m^3}}$",size=27)
plt.legend()
plt.savefig("rho120sedov.pdf")
plt.close()

plt.figure(figsize=(20,20))
plt.plot(x,rho2)
plt.text(20.0,1.1,"t="+str(dt2)+"s",va='center',ha='left',size=20)
plt.xlim(0,128)
plt.xlabel(r"$\mathrm{x}$ $(m)$",size=20)
plt.ylabel(r"$\mathrm{\rho}$ $\mathrm{\frac{Kg}{m^3}}$",size=27)
plt.legend()
plt.savefig("rho60sedov.pdf")
plt.close()

plt.figure(figsize=(20,20))
plt.plot(x,rho1)
plt.text(20.0,1.1,"t="+str(dt1)+"s",va='center',ha='left',size=20)
plt.xlim(0,128)
plt.xlabel(r"$\mathrm{x}$ $(m)$",size=20)
plt.ylabel(r"$\mathrm{\rho}$ $\mathrm{\frac{Kg}{m^3}}$",size=27)
plt.legend()
plt.savefig("rho10sedov.pdf")
plt.close()


