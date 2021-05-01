# Code for solving the Gross-Pitaevskii equation in 2D and 3D using parallel computing employing multiprocessing. We use the explicit 
# Euler forward method for time evolution. The array is divided along z-axis for the purpose of parallelising the code.
# Number of processors to be used has to be provided as command line argument. E.g. use 'python GP_mp 2' for using 2 processors in multiprocessing.

import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.integrate import trapz
import multiprocessing
import SharedArray

# To make 2D colormap plot of the solution (|\Psi|^2 as a function of x and z) at given time t. The figure is saved
# at the location given in variable filename. For 3D code, we plot at a particular value of y coordinate.
# The variables x and z contain points along x and z, and prob contains |\Psi|^2.
def plot2D(x, z, prob, t, i):
	plt.close()
	fig = plt.figure()
	plt.axis('square')
	ax = plt.gca()
	probmap = ax.pcolormesh(x, z, np.transpose(prob), cmap = cm.jet, shading = 'auto', vmin = 0.0, vmax = 1.2)
	cbar = plt.colorbar(probmap, orientation='vertical')
	cbar.set_label(r'$|\psi^2|$', fontsize = 14, rotation = 0, labelpad = 20)
	ax.set_title(r't = %.9f' %(t))
	ax.set_xlabel(r'$x$')
	ax.set_ylabel(r'$z$')
	ax.set_xticks([-1, 0, 1])
	ax.set_yticks([-1, 0, 1])
	ax.set_xlim([x[0], x[-1]])
	ax.set_ylim([z[0], z[-1]])
	filename = 'path/fig_%03d.png' %(i)
	plt.savefig(filename)
	
	return

# Set initial condition for Psi. The initial condition must satisfy the boundary conditions,
# which is \Psi = 0 at all the boundaries and at all times. x, y and z variables must be output
# of np.meshgrid function (if normal 1D arrays are used, the output value of Psi will also be a 1D array).
def init_cond(x, y, z):
	global ny
	# 2D
	if ny == 1:
		init_Psi = np.sin(np.pi*x)*np.sin(np.pi*z)*np.ones_like(y)
	# 3D
	else:
		init_Psi = np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(2*np.pi*z)

	return -1j*1j*init_Psi

# Potential. x, y and z variables must be output of np.meshgrid function
# (if normal 1D arrays are used, the output value of Psi will also be a 1D array).
def V(x, y, z):
	global ny
	# 2D
	if ny == 1:
		potential = 400*(0.5*x**2 + 0.5*z**2)
	# 3D
	else:
		potential = 400*(0.5*x**2 + 0.5*y**2 + 0.5*z**2)

	return potential

# This function computes the evolution of wavefunction Psi for 2D and again stores it in the variable Psi.
# The array Psi is shared across all the processes using SharedArray.
def H_2D(proc_id):
	# The start and end index along z axis that processor with id = proc_id has to compute.
	start = proc_id*arr_len + 1
	end = (proc_id + 1)*arr_len + 1

	# var0 computes the part of the Hamiltonian containing potential and non-linear term.
	var0 = (potential[1:nx-1, :, start:end] + N*np.abs(Psi[1:nx-1, :, start:end])**2)*Psi[1:nx-1, :, start:end]

	# var1 computes part of the potential containg second order derivatie along x.
	var1 = (Psi[2:nx, :, start:end] - 2*Psi[1:nx-1, :, start:end] + Psi[0:nx-2, :, start:end])/dx**2

	# var3 computes part of the potential containg second order derivatie along z.
	var3 = (Psi[1:nx-1, :, start+1:end+1] - 2*Psi[1:nx-1, :, start:end] + Psi[1:nx-1, :, start-1:end-1])/dz**2

	# Computing time evolution using Euler forward.
	Psi[1:nx-1, :, start:end] = Psi[1:nx-1, :, start:end] + dt*1j*(var1 + var3 - var0)

	return

# This function computes the evolution of wavefunction Psi for 3D and again stores it in the variable Psi.
# The array Psi is shared across all the processes using SharedArray.
def H_3D(proc_id):
	# The start and end index along z axis that processor with id = proc_id has to compute.
	start = proc_id*arr_len + 1
	end = (proc_id + 1)*arr_len + 1

	# var0 computes the part of the Hamiltonian containing potential and non-linear term.
	var0 = (potential[1:nx-1, 1:ny-1, start:end] + N*np.abs(Psi[1:nx-1, 1:ny-1, start:end])**2)*Psi[1:nx-1, 1:ny-1, start:end]

	# var1 computes part of the potential containg second order derivatie along x.
	var1 = (Psi[2:nx, 1:ny-1, start:end] - 2*Psi[1:nx-1, 1:ny-1, start:end] + Psi[0:nx-2, 1:ny-1, start:end])/dx**2

	# var2 computes part of the potential containg second order derivatie along y.
	var2 = (Psi[1:nx-1, 2:ny, start:end] - 2*Psi[1:nx-1, 1:ny-1, start:end] + Psi[1:nx-1, 0:ny-2, start:end])/dy**2

	# var3 computes part of the potential containg second order derivatie along z.
	var3 = (Psi[1:nx-1, 1:ny-1, start+1:end+1] - 2*Psi[1:nx-1, 1:ny-1, start:end] + Psi[1:nx-1, 1:ny-1, start-1:end-1])/dz**2

	# Computing time evolution using Euler forward.
	Psi[1:nx-1, 1:ny-1, start:end] = Psi[1:nx-1, 1:ny-1, start:end] + dt*1j*(var1 + var2 + var3 - var0)

	return

# Computes the evolution of Psi with respect to t.
def evolve(x, y, z, t, skip, nprocs):
	# Inputs array containing a list of process ids.
	inputs = list(range(nprocs))
	pool = multiprocessing.Pool(nprocs)

	# Differentiaing between 2D and 3D code.
	if ny == 1:
		func = H_2D
	else:
		func = H_3D

	# Loop for computing time evolution
	start = time.time()
	for k in range(1, m):
		# Plotting and computing integral at specific time steps. Uncomment this section if you want to output plots and integral values.
		'''
		if k % skip == 1:
			if ny == 1:
				plot2D(x, z, np.abs(Psi[:, 0, :])**2, t[k], int(k/skip))
				integral[int(k/skip)] = trapz(trapz(np.abs(Psi[:, 0, :])**2, z), x)
			else:
				plot2D(x, z, np.abs(Psi[:, ny//4, :])**2, t[k], int(k/skip))
				integral[int(k/skip)] = trapz(trapz(trapz(np.abs(Psi)**2, z), y), x)
		'''

		pool.map(func, inputs)

	stop = time.time()

	# Printing the time taken.
	print('Total Time:', stop - start)
	return

# Number of processes to be used. Has to be given as command line argument
# For example, use 'python GP_mp.py 2' to use nprocs = 2.
nprocs = int(sys.argv[1])

# Time step (use a small time step to avoid instability)
dt = 1e-9

# Length of x, y, z, t arrays
# 2 added to nx, ny and nz is for the first and last z index at which Psi = 0 due to boundary condition.
nx = 1024 + 2
ny = 1
nz = 1024 + 2
m = 10001
arr_len = int((nz-2)/nprocs)

# Parameter to be multiplied with the non-linear term
N = 10.0

# x, y, z, t arrays
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
z = np.linspace(-1, 1, nz)
t = np.linspace(0, (m-1)*dt, m)

# Difference between two successive x, y, z, and t values, respectively.
dx = x[1] - x[0]
dz = z[1] - z[0]
if ny != 1:
	dy = y[1] - y[0]
dt = t[1] - t[0]

# Meshgrid
xv, yv, zv = np.meshgrid(x, y, z, indexing = 'ij')

# Potential
potential = V(xv, yv, zv)

if __name__ == "__main__":
    # Computes the interval at which to plot colormap and compute integral.
    framerate = 24
    vidlen = 1
    no_of_frames = framerate*vidlen
    skip = int(np.floor((m - 1)/no_of_frames))

    # Array to compute integral of |\Psi|^2
    integral = np.zeros(int(np.floor(m/skip) + 1))
    indices = skip*np.linspace(0, len(integral) - 1, len(integral), dtype = np.int)
    t_arr = t[indices]

    # SharedArray to store Psi.
    a = SharedArray.create("shm://Psi", (nx, ny, nz), dtype = np.complex128)

# Attaching the SharedArray to the variable Psi for each process.
Psi = SharedArray.attach("shm://Psi")
# Initialising Psi
Psi[:] = init_cond(xv, yv, zv)

if __name__ == "__main__":
	# Evolution
    evolve(x, y, z, t, skip, nprocs)

    # Deleting SharedArray
    SharedArray.delete("shm://Psi")

    # Plotting integral vs t. Uncomment this section if you are computing the integral inside evolve function.
    '''
    plt.close()
    plt.plot(t_arr, integral)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\int|\psi|^2\ \mathregular{d}^2x$')
    plt.savefig('figures/int_par.png', dpi = 300)
    '''



