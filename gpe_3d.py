import time
import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

def mat(v):
    if v=='x':
        (A,B,C,N)=(Ax,Bx,Cx,Nx)
    elif v=='y':
        (A,B,C,N)=(Ay,By,Cy,Ny)
    else:
        (A,B,C,N)=(Az,Bz,Cz,Nz)
    M=np.zeros((N-2,N),dtype=complex)
    for i in range(0,N-2):
        M[i,i]=A
        M[i,i+1]=B
        M[i,i+2]=C
    return(M[:,1:-1])


def H1(temppsi):
    r=np.zeros((Nz-2,Ny-2,Nx-2),dtype=complex)
    for i in range(0,Nx-2):
        r[:,:,i]=-Ax*(temppsi[1:-1,1:-1,i+2]-2*temppsi[1:-1,1:-1,i+1]+temppsi[1:-1,1:-1,i])+temppsi[1:-1,1:-1,i+1]
    return(np.transpose(solve(mat('x'),np.transpose(r,(2,1,0))),(2,1,0)))

def H2(temppsi):
    r=np.zeros((Nz-2,Ny-2,Nx-2),dtype=complex)
    for i in range(0,Ny-2):
        r[:,i,:]=-Ay*(temppsi[1:-1,i+2,1:-1]-2*temppsi[1:-1,i+1,1:-1]+temppsi[1:-1,i,1:-1])+temppsi[1:-1,i+1,1:-1]
    return(np.transpose(solve(mat('y'),np.transpose(r,(1,0,2))),(1,0,2)))

def H3(temppsi):
    r=np.zeros((Nz-2,Ny-2,Nx-2),dtype=complex)
    for i in range(0,Nz-2):
        r[i,:,:]=-Az*(temppsi[i+2,1:-1,1:-1]-2*temppsi[i+1,1:-1,1:-1]+temppsi[i,1:-1,1:-1])+temppsi[i+1,1:-1,1:-1]
    return(solve(mat('z'),r))

def initfunc():
    return 2/l*np.sin(np.pi*xv/l)*np.sin(np.pi*yv/l)#*np.sin(np.pi*zv/l);

c=complex(0,1)#complex number
l=6 #length of infinite square well
Tt=10#total time in seconds
n=1001#total number of timing points including t=0
Nx=32#total no of grid points along x
Ny=32 #total no of grid points along y
Nz=32#total no of grid points along z
x=np.linspace(0,l,Nx)
y=np.linspace(0,l,Ny)
z=np.linspace(0,l,Nz)
t=np.linspace(0,Tt,n)
dx=l/(Nx-1)
dy=l/(Ny-1)
dz=l/(Nz-1)
dt=Tt/(n-1)

yv,zv,xv=np.meshgrid(y,z,x)
Ax=complex(0,-dt/(2*dx**2))
Bx=complex(1,dt/dx**2)
Cx=Ax

Ay=complex(0,-dt/(2*dy**2))
By=complex(1,dt/dy**2)
Cy=Ay

Az=complex(0,-dt/(2*dz**2))
Bz=complex(1,dt/dz**2)
Cz=Az


psi=np.zeros((n,Nz,Ny,Nx),dtype=complex)  #wavefunction each raw contain data for each point in space and for constant time
temppsi=np.zeros((Nz,Ny,Nx),dtype=complex)  # contain n+1/2 value of the wavefunction

psi[0,:,:,:]=initfunc()   #initial function
for i in range(n-1):
    temppsi=(1-dt*abs(psi[i,:,:,:])**2*c)*psi[i,:,:,:]
    psi[i+1,1:-1,1:-1,1:-1]=H1(temppsi)
    psi[i+1,1:-1,1:-1,1:-1]=H2(psi[i+1,:,:,:])
    psi[i+1,1:-1,1:-1,1:-1]=H3(psi[i+1,:,:,:])
