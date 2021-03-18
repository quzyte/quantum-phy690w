import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt


def crank_nik(temppsi):
    r=np.zeros(Nx-2,dtype=complex)
    M=np.zeros((Nx-2,Nx),dtype=complex)
    for i in range(0,Nx-2):
        r[i]=-A*(temppsi[i+2]-2*temppsi[i+1]+temppsi[i])+temppsi[i+1]
        M[i,i]=A
        M[i,i+1]=B
        M[i,i+2]=C
    M=M[:,1:-1]
    return(solve(M, r))


def initfunc():
    return np.sqrt(2/l)*np.sin(3*np.pi*x/l);

c=complex(0,1)#complex number
l=6    #length of infinite square well
Tt=1   #total time in seconds
n=101  #total number of timing points including t=0
Nx=101 #total no of grid points

x=np.linspace(0,l,Nx)
t=np.linspace(0,Tt,n)
dx=l/(Nx-1)
dt=Tt/(n-1)

A=complex(0,-dt/(2*dx**2))
B=complex(1,dt/dx**2)
C=A
psi=np.zeros((n,Nx),dtype=complex)  #wavefunction each raw contain data for each point in space and for constant time
temppsi=np.zeros(Nx,dtype=complex)  # contain n+1/2 th value of the wavefunction

psi[0,:]=initfunc()   #initial function


for i in range(n-2):
    temppsi=(1-dt*abs(psi[i,:])**2*c)*psi[i,:]
    psi[i+1,1:-1]=crank_nik(temppsi)
for i in range(0,n,10):
    plt.plot(x,np.real(psi[i,:]))
plt.show()
