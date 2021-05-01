
# Code for solving the Gross-Pitaevskii equation in 2D and 3D using parallel computing employing MPI. We use the explicit 
# Euler forward method for time evolution. The array is divided along z-axis for the purpose of parallelising the code.

# The size of x, y and z arrays used are: 2^n + 2. Suppose that we choose n = 10, i.e. size = 1024 + 2, and we are using 4 processors.
# Then, array size along z for each processor will then be 256 + 2 = 258, including 2 dummy indices. Dummy indices are the first and last index
# along z which are used to compute the differentiation along z, but value of Psi at next time step is not computed at these dummy indices by the 
# processor. The value of Psi at these dummy indices is updated after communication between processes. For first processor, value of Psi at first z
# index is 0 (irrespective of x and y index) because of boundary condition. So, we can treat this as dummy index as well. Similarly, for last
# processor, Psi = 0 at last z index.

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from scipy.integrate import trapz, simps
from scipy.linalg import solve_banded
from mpi4py import MPI

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

# This function computes the -i H \Psi for 2D problem. Here, H is the Hamiltonian operator, and i is sqrt(-1).
# The output of this function is related to time derivative of \Psi by d\Psi/dt = output
def H_2D(Psi, potential, dx, dy, dz, nx, ny, nz, N):
	# var0 computes the part of the Hamiltonian containing potential and non-linear term.
	var0 = ((potential + N*np.abs(Psi)**2)*Psi)[1:nx-1, 0, 1:nz-1]

	# var1 computes part of the potential containg second order derivatie along x.
	var1 = (Psi[2:nx, 0, 1:nz-1] - 2*Psi[1:nx-1, 0, 1:nz-1] + Psi[0:nx-2, 0, 1:nz-1])/dx**2

	# var3 computes part of the potential containg second order derivatie along z.
	var3 = (Psi[1:nx-1, 0, 2:nz] - 2*Psi[1:nx-1, 0, 1:nz-1] + Psi[1:nx-1, 0, 0:nz-2])/dz**2

	return 1j*(var1 + var3 - var0)

# This function computes the -i H \Psi for 2D problem. Here, H is the Hamiltonian operator, and i is sqrt(-1).
# The output of this function is related to time derivative of \Psi by d\Psi/dt = output
def H_3D(Psi, potential, dx, dy, dz, nx, ny, nz, N):
	# var0 computes the part of the Hamiltonian containing potential and non-linear term.
	var0 = ((potential + N*np.abs(Psi)**2)*Psi)[1:nx-1, 1:ny-1, 1:nz-1]

	# var1 computes part of the potential containg second order derivatie along x.
	var1 = (Psi[2:nx, 1:ny-1, 1:nz-1] - 2*Psi[1:nx-1, 1:ny-1, 1:nz-1] + Psi[0:nx-2, 1:ny-1, 1:nz-1])/dx**2

	# var2 computes part of the potential containg second order derivatie along y.
	var2 = (Psi[1:nx-1, 2:ny, 1:nz-1] - 2*Psi[1:nx-1, 1:ny-1, 1:nz-1] + Psi[1:nx-1, 0:ny-2, 1:nz-1])/dy**2

	# var3 computes part of the potential containg second order derivatie along z.
	var3 = (Psi[1:nx-1, 1:ny-1, 2:nz] - 2*Psi[1:nx-1, 1:ny-1, 1:nz-1] + Psi[1:nx-1, 1:ny-1, 0:nz-2])/dz**2

	return 1j*(var1 + var2 + var3 - var0)

# Computes the evolution of Psi with respect to t.
def evolve(nx, ny, nz, m, dt, N, nprocs, rank, skip):
	# x, y, z, t arrays. We're using x, y, z arrays from -1 to +1 because we are using quadratic potential,
	x = np.linspace(-1, 1, nx)
	y = np.linspace(-1, 1, ny)
	z = np.linspace(-1, 1, nz)
	t = np.linspace(0, (m-1)*dt, m)

	# Difference between two successive x, y, z, and t values, respectively.
	dx = x[1] - x[0]
	if ny == 1:
		# For 2D code, array along y does not matter.
		dy = 0
	else:
		dy = y[1] - y[0]
	dz = z[1] - z[0]
	dt = t[1] - t[0]

	# Length of the z array per processor.
	arr_len = int((nz-2)/nprocs)

	# z array for each processor.
	z_part = z[rank*arr_len:(rank+1)*arr_len + 2]

	# No. of elements in z array for each processor.
	nz = arr_len + 2

	# Meshgrid
	xv, yv, zv = np.meshgrid(x, y, z_part, indexing = 'ij', sparse = True)

	# Initialising Psi and potential
	Psi = init_cond(xv, yv, zv)
	potential = V(xv, yv, zv)

	# Arrays to be used for communication among processors.
	send_buf = np.zeros_like(Psi[:, :, 0], dtype = np.complex128)
	recv_buf = np.zeros_like(Psi[:, :, 0], dtype = np.complex128)

	# Array we will use to get data from each processor to rank 0 for the purpose of plotting and computing integral.
	Psi_temp = np.zeros_like(Psi[:, :, 2:], dtype = np.complex128)

	# Differentiaing between 2D and 3D code.
	if ny == 1:
		y_range = 0
		func = H_2D
	else:
		y_range = np.arange(1, ny-1)
		func = H_3D

	# Computational time and Communication time
	comp_time = 0.0
	comm_time = 0.0

	# For communication, we use the following strategy. First, we pair the processors as follows: 0-1, 2-3, 4-5, 6-7, ...
	# We use Sendrecv between these pairs of processors. So, an even ranked processor first communicates with the next processor
	# (i.e. 0 communicates with 1, 2 with 3 and so on), while an odd ranked processor communicates with previous one.
	# Then, we pair processors as: 1-2, 3-4, 5-6, ...  So, even ranked communicates with previous one and odd rank communicates
	# with next one. This information in encompassed in the variables first_rank and second_rank.
	# The variables first_send_index and second_send_index contains the z index to be sent during first and second set of
	# communication. Similarly, first_recv_index and second_recv_index are the z indices at which recv_bufs will go.
	if rank % 2 == 0:
		first_rank = rank + 1
		first_send_index = nz-2
		first_recv_index = nz-1

		second_rank = rank - 1
		second_send_index = 1
		second_recv_index = 0

	else:
		first_rank = rank - 1
		first_send_index = 1
		first_recv_index = 0

		second_rank = rank + 1
		second_send_index = nz-2
		second_recv_index = nz-1

	# For first processor and last processor, second_rank should be -1, so that NULL communication occurs. This would prevent
	# MPI from raising error.
	if rank == nprocs-1:
		second_rank = -1
	if nprocs == 1:
		first_rank = -1

	# Loop for computing time evolution
	for k in range(1, m):
		# Plotting and computing integral at specific time steps. Uncomment this section if you want to output plots and integral values.
		# Here, first each processor will send their Psi arrays to rank 0, where the complete Psi array will be built. The integral and
		# plotting will then be done by rank 0.
		'''
		if k % skip == 1:
			if rank == 0:
				Psi_tot = Psi
				for i in range(1, nprocs):
					comm.Recv([Psi_temp, MPI.COMPLEX], source = i, tag = i)
					Psi_tot = np.concatenate((Psi_tot, Psi_temp), axis = 2)
				
				# 2D
				if ny == 1:
					plot2D(x, z, np.abs(Psi_tot[:, 0, :])**2, t[k], int(k/skip))
					integral[int(k/skip)] = trapz(trapz(np.abs(Psi_tot[:, 0, :])**2, z), x)
				# 3D
				else:
					plot2D(x, z, np.abs(Psi_tot[:, ny//4, :])**2, t[k], int(k/skip))
					integral[int(k/skip)] = trapz(trapz(trapz(np.abs(Psi_tot)**2, z), y), x)

			else:
				comm.Send([np.ascontiguousarray(Psi[:, :, 2:]), MPI.COMPLEX], dest = 0, tag = rank)
		'''

		# Computation
		start_comp = time.time()
		# Euler forward method
		Psi[1:nx-1, y_range, 1:nz-1] = Psi[1:nx-1, y_range, 1:nz-1] + dt*func(Psi, potential, dx, dy, dz, nx, ny, nz, N)
		
		end_comp = time.time()
		comp_time = comp_time + end_comp - start_comp

		# Communication
		start_comm = time.time()
		# First communication set
		send_buf = np.ascontiguousarray(Psi[:, :, first_send_index])
		recv_buf = 0*recv_buf
		comm.Sendrecv([send_buf, MPI.COMPLEX], dest = first_rank, sendtag = 100, recvbuf = [recv_buf, MPI.COMPLEX], source = first_rank, recvtag = 100)
		Psi[:, :, first_recv_index] = recv_buf
		# Second communication set
		send_buf = np.ascontiguousarray(Psi[:, :, second_send_index])
		recv_buf = 0*recv_buf
		comm.Sendrecv([send_buf, MPI.COMPLEX], dest = second_rank, sendtag = 200, recvbuf = [recv_buf, MPI.COMPLEX], source = second_rank, recvtag = 200)
		Psi[:, :, second_recv_index] = recv_buf

		end_comm = time.time()
		comm_time = comm_time + end_comm - start_comm

	# Printing the time taken.
	if rank == 0:
		print('Communication Time:', comm_time)
		print('Computation Time:', comp_time)
		print('Total Time:', comm_time + comp_time)

	return


# Variables required for MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

# Time step (use a small time step to avoid instability)
dt = 1e-9

# Length of x, y, z, t arrays
nx = 1024 + 2
ny = 1
nz = 1024 + 2
m = 1001

# Parameter to be multiplied with the non-linear term
N = 10.0

# Computes the interval at which to plot colormap and compute integral.
framerate = 24
vidlen = 1
no_of_frames = framerate*vidlen
skip = int(np.floor((m - 1)/no_of_frames))

# Array to compute integral of |\Psi|^2
integral = np.zeros(int(np.floor(m/skip) + 1))
indices = skip*np.linspace(0, len(integral) - 1, len(integral), dtype = np.int)
t = np.linspace(0, (m-1)*dt, m)
t_arr = t[indices]

# Function which computes evolution
evolve(nx, ny, nz, m, dt, N, nprocs, rank, skip)

# Plotting integral vs t. Uncomment this section if you are computing the integral inside evolve function.
'''
if rank == 0:
	plt.close()
	plt.plot(t_arr, integral)
	plt.xlabel(r'$t$')
	plt.ylabel(r'$\int|\psi|^2\ \mathregular{d}^2x$')
	plt.savefig('figures/int_par.png', dpi = 300)
'''

