from matplotlib import pyplot as plt
import numpy as np


#HEMEN UNITATE NATURALAK ERABILTZEN DITU BAINA GUK ALDATU DEZAKEGU
hbar=1
m=1
omega=1
k=2

#N TARTE KOPURUA IZANGO DA ETA A GURE TARTEAREN LUZEERA 
N = 1000
a = 20.0

#HAU FORTRANEN "TARTA" DA
x = np.linspace(-a/2.,a/2.,N)

#ETA HAU GURE PROGRAMAKO H ERE BAI BAINA BESTE ERAN DEFINITUTA
h = x[1]-x[0] # Should be equal to 2*np.pi/(N-1)

#POTENTZIALA FINKATUKO BEHAR DUGU
V2 = .5*k*x*x
V1=np.full(N,0)
esk=100
E=100
a=E/esk
E_bek=np.full(N,E)

for i in range(len(x)):
    S=0
    for j in range(esk):
        muga=.5*k*x[i]*x[i]/a+0.5
        if muga>j:
            V1[i]=V1[i]+a
        elif muga<j:
            V1[i]=V1[i]

#AQUI ESTA EL MEOLLO :) D HACE LA MATRIZ DE LAS DERIVADAS Y EN H LAS SUMA Y DESPUES SOLO QEUDA DIAGONALIZAR
D= 1./(h*h)*(np.diag(np.ones(N-1),-1) -2* np.diag(np.ones(N),0) + np.diag(np.ones(N-1),1))
H1 = -(hbar*hbar)/(2.0*m)*D + np.diag(V1)
H2 = -(hbar*hbar)/(2.0*m)*D + np.diag(V2)

#ETA HEMEN LA MAGIA DE DIAGONALIZAR
En1,psiT1 = np.linalg.eigh(H1) # This computes the eigen values and eigenvectors
psi1 = np.transpose(psiT1)

En2,psiT2=np.linalg.eigh(H2)
psi2 = np.transpose(psiT2)



#HEMEN GRAFIKOA EGINGO DUGU


fig, ax1 = plt.subplots()

ax1.set_xlabel('$x$')
ax1.set_ylabel('$\psi(x)$')

ax2 = ax1.twinx()

ax2.set_ylabel('$V(x)$')
plt.title('Harmonic Oscillator')



ax1.grid(True, linestyle='-.')
plt.xlim((-10.,10.))
for i in range(0,3):
    if psi1[i][N//8] < 0:
        ax1.plot(x,-psi1[i]/np.sqrt(h),label="$E_{}$={:3.1f}".format(i,En1[i]))
    else:
        ax1.plot(x,psi1[i]/np.sqrt(h),label="$E_{}$={:3.1f}".format(i,En1[i]))
    if psi2[i][N//8] < 0:
        ax1.plot(x,-psi2[i]/np.sqrt(h),label="$E_{}$={:3.1f}".format(i,En2[i]),linestyle='dashed',color="grey",linewidth=0.5)
    else:
        ax1.plot(x,psi2[i]/np.sqrt(h),label="$E_{}$={:3.1f}".format(i,En2[i]),linestyle='dashed',color="grey",linewidth=0.5)

ax2.plot(x,V1,color="blue",label="$V_1(x)$")
ax2.plot(x,V2,color="grey",label="$V_2(x)$",linewidth=0.5,linestyle='dashed')
ax2.plot(x,E_bek,color="red",label="$E$")

plt.title("Solution to harmonic oscillator")
plt.legend()
plt.savefig("Harmonic_Oscillator_WaveFunctions.svg")
plt.show()


