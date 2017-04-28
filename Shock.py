import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
u=np.loadtxt("ushock.dat")
p=np.loadtxt("pshock.dat")
rho=np.loadtxt("rhoshock.dat")
t=np.loadtxt("tshock.dat")
x=np.linspace(0,1,len(u))
x43=np.linspace(0,0.2498,len(u))
x32p=np.linspace(0.2498,0.498,498-250)
x22=np.linspace(0.498,0.70,len(u))
x21=np.linspace(0.70,0.9,len(u))
xf=np.linspace(0.9,1.0,len(u))
xi=np.linspace(0.498,0.9,len(u))




a4=-1.183215957
pz1=[]
pz2=[]
pz3=[]
pz4=[]
uz1=[]
uz2=[]
uz4=[]
rhoz1=[]
rhoz2=[]
rhoz2p=[]
rhoz3=[]
rhoz4=[]
for i in range(len(u)):
    pz4.append(1.0)
    pz2.append(0.284)
    pz1.append(0.1)
    uz1.append(0.0)
    uz2.append((a4/1.4)*(pz2[0]-1.0)*np.sqrt((2.0*1.4/2.4)/(pz2[0]+(0.4/2.4))))
    uz4.append(0.0)
    rhoz2.append(0.2044)
    rhoz2p.append(0.407)
    rhoz4.append(1.0)
    rhoz1.append(0.1)
    if(i>=250 and i<498):
         rhoz3.append(p[i]**(1.0/1.4))
         pz3.append(rho[i]**1.4)
   
def uz3(x):
    return (2.0/2.4)*(a4+(x/t))



plt.figure(figsize=(20,20))
plt.plot(x,u,"k",label="Lax-wendroff")
plt.plot(x32p,uz3(x32p),"--",label="Teorico")
plt.plot(xf,uz1,"--")
plt.plot(xi,uz2,"--")
plt.plot(x43,uz4,"--")
plt.ylim(-0.1,1.1)
plt.xlabel(r"$\mathrm{X(m)}$",size=20)
plt.ylabel(r"$\mathrm{v}$ $(\mathrm{\sqrt{101000}\frac{m}{s}}) $",size=20)
plt.legend()
plt.savefig("u.pdf")
plt.close()
plt.figure(figsize=(20,20))
plt.plot(x,p,"k",label="Lax-wendroff")
plt.plot(xf,pz1,"--",label="Teorico")
plt.plot(xi,pz2,"--")
plt.plot(x32p,pz3,"--")
plt.plot(x43,pz4,"--")
plt.ylim(0.0,1.1)
plt.xlabel(r"$\mathrm{X(m)}$",size=20)
plt.ylabel(r"$\mathrm{P}$ $(atm)$",size=20)
plt.legend()
plt.savefig("p.pdf")
plt.close()
plt.figure(figsize=(20,20))
plt.plot(x,rho,"k",label="Lax-wendroff")
plt.plot(xf,rhoz1,"--",label="Teorico")
plt.plot(x21,rhoz2,"--")
plt.plot(x22,rhoz2p,"--")
plt.plot(x32p,rhoz3,"--")
plt.plot(x43,rhoz4,"--")
plt.ylim(0.0,1.1)
plt.xlabel(r"$\mathrm{X(m)}$",size=20)
plt.ylabel(r"$\mathrm{\rho}$ $(\frac{Kg}{m^3})$",size=27)
plt.legend()
plt.savefig("rho.pdf")
plt.close()



