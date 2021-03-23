import time
import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.cm as cm
from pylab import rcParams
rcParams['figure.figsize'] = 5, 5



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
    r=np.zeros((Ny-2,Nx-2),dtype=complex)
    for i in range(0,Nx-2):
        r[:,i]=-Ax*(temppsi[1:-1,i+2]-2*temppsi[1:-1,i+1]+temppsi[1:-1,i])+temppsi[1:-1,i+1]
    return(np.transpose(solve(mat('x'),np.transpose(r))))

def H2(temppsi):
    r=np.zeros((Ny-2,Nx-2),dtype=complex)
    for i in range(0,Ny-2):
        r[i,:]=-Ay*(temppsi[i+2,1:-1]-2*temppsi[i+1,1:-1]+temppsi[i,1:-1])+temppsi[i+1,1:-1]
    return(solve(mat('y'),r))
'''
#def H3(temppsi):
#    sc
'''
def initfunc():
    return 2/l*np.sin(np.pi*xv/l)*np.sin(np.pi*yv/l)#*np.sin(np.pi*zv/l);

c=complex(0,1)#complex number
l=6 #length of infinite square well
Tt=10#total time in seconds
n=1001#total number of timing points including t=0
Nx=105 #total no of grid points along x
Ny=101 #total no of grid points along y
Nz=101 #total no of grid points along z
x=np.linspace(0,l,Nx)
y=np.linspace(0,l,Ny)
z=np.linspace(0,l,Nz)
t=np.linspace(0,Tt,n)
dx=l/(Nx-1)
dy=l/(Ny-1)
dz=l/(Nz-1)
dt=Tt/(n-1)

xv,yv=np.meshgrid(x,y)
Ax=complex(0,-dt/(2*dx**2))
Bx=complex(1,dt/dx**2)
Cx=Ax

Ay=complex(0,-dt/(2*dy**2))
By=complex(1,dt/dy**2)
Cy=Ay

Az=complex(0,-dt/(2*dz**2))
Bz=complex(1,dt/dz**2)
Cz=Az


psi=np.zeros((n,Ny,Nx),dtype=complex)  #wavefunction each raw contain data for each point in space and for constant time
temppsi=np.zeros((Ny,Nx),dtype=complex)  # contain n+1/2 value of the wavefunction

psi[0,:,:]=initfunc()   #initial function
for i in range(n-1):
    temppsi=(1-dt*abs(psi[i,:,:])**2*c)*psi[i,:,:]
    psi[i+1,1:-1,1:-1]=H1(temppsi)
    psi[i+1,1:-1,1:-1]=H2(psi[i+1,:,:])
phi=np.real(psi)
'''
fig = plt.figure()
axes = fig.gca(projection='3d')
p = axes.plot_surface(xv,yv,phi[0,:,:], rstride=4,
cstride=4, cmap=cm.RdBu, linewidth=0, antialiased=False)
fig.colorbar(p, shrink=0.5)
axes.set_xlabel('$x$',fontsize=15)
axes.set_ylabel('$y$',fontsize=15)
axes.set_zlabel('$z$',fontsize=15)
plt.tight_layout();
plt.show()

'''
fps=100
def update_plot(frame_number,phi, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(xv, yv,phi[frame_number,:,:], cmap="magma")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plot = [ax.plot_surface(xv, yv,phi[0,:,:], color='0.75', rstride=1, cstride=1)]
ax.set_zlim(0,.3)
ani = animation.FuncAnimation(fig, update_plot, n, fargs=(phi, plot), interval=10)
fn = 'plot_surface_animation_funcanimation'

ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)
#ani.save(fn+'.gif',writer='imagemagick',fps=fps)
