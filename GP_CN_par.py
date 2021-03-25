import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from scipy.integrate import trapz
from scipy.linalg import solve_banded
from mpi4py import MPI

# Set initial condition for Psi
def init_cond(x, y, z):
	init_Psi = np.zeros([len(x), len(y), len(z)], dtype = np.complex128)

	# 2D
	if ny == 1:
		init_Psi = 2*np.sin(np.pi*x)*np.sin(2*np.pi*z)*np.ones_like(y)

	# 3D
	else:
		init_Psi = np.sin(np.pi*x)*np.sin(2*np.pi*y)*np.sin(3*np.pi*z)

	return init_Psi

# 2D Potential
def V_2D(x, y, z, nx, ny, nz):
	ret_val = np.zeros([nx, ny, nz], dtype = np.complex128)
	ret_val = 40*(0.5*x**2 + 0.5*z**2)

	return ret_val

# 3D Potential
def V_3D(x, y, z, nx, ny, nz):
	ret_val = np.zeros([nx, ny, nz], dtype = np.complex128)
	ret_val = 0.5*x**2 + 0.5*y**2 + 0.5*z**2

	return ret_val

# Returns the solution to 2D equation with H_0 as Hamiltonian. Does not involve any spatial derivative
def H0_2D(Psi, x, y, z, nx, ny, nz, dt, N):
	potential = V_2D(x, y, z, nx, ny, nz) + N*np.abs(Psi)**2
	return np.exp(-1j*dt*potential)*Psi

# Returns the solution to 3D equation with H_0 as Hamiltonian. Does not involve any spatial derivative
def H0_3D(Psi, x, y, z, nx, ny, nz, dt, N):
	potential = V_3D(x, y, z, nx, ny, nz) + N*np.abs(Psi)**2
	return np.exp(-1j*dt*potential)*Psi

# For Functions H1, H2, H3:
# We have to solve matrix equation Ax = b for x, where A is a tridiagonal matrix.
# However, here, b is not a simple column matrix. Instead, it's a 3D matrix. 
# We have to solve the equation Ax = b for each column in b.
# For eample, when we're solving the equation with H_1 (that is x derivative), we need to solve Ax = b[:, i, j] for all i and j
# We can do this using for loops. But, for loop will slow down the code.
# So, we use a scipy function solve_banded(). solve_banded() allows us to pass b as a 3D matrix
# It will return the solution x, also a 3D matrix of ths same shape a b.
# However, solve_banded solves the equation Ax = b along axis 0 only.
# That is, it will return the solution A*x[:, i, j] = b[:, i, j] for all i and j.
# This would cause a problem when we'ra solving equation with H_2 or H_3 (that is, derivative with respect to y or z).
# To overcome this problem, we first transpose b in such a way that it is appropriate to pass b directly to solve_banded.
# Later, we would again need to transpose the solution.

# Returns the solution to equation with H_1 as Hamiltonian. Involves derivative w.r.t. x
# Banded_matrix is in the form requrired by the function solve_banded(). See scipy documentation.
def H1(Psi, banded_matrix, mu, nx, ny, nz):
	b = mu*Psi[2:nx, :, :] + (1-2*mu)*Psi[1:nx-1, :, :] + mu*Psi[0:nx-2, :, :]

	Psi_new = np.zeros([nx, ny, nz], dtype = np.complex128)
	Psi_new[1:nx-1, :, :] = solve_banded((1, 1), banded_matrix, b)
	
	return Psi_new

# Returns the solution to equation with H_1 as Hamiltonian. Involves derivative w.r.t. y
def H2(Psi, banded_matrix, mu, nx, ny, nz):
	b = mu*Psi[:, 2:ny, :] + (1-2*mu)*Psi[:, 1:ny-1, :] + mu*Psi[:, 0:ny-2, :]
	b = np.transpose(b, axes = [1, 0, 2])

	Psi_new = np.zeros([ny, nx, nz], dtype = np.complex128)
	Psi_new[1:ny-1, :, :] = solve_banded((1, 1), banded_matrix, b)
	Psi_new = np.transpose(Psi_new, axes = [1, 0, 2])
	
	return Psi_new

# Returns the solution to equation with H_1 as Hamiltonian. Involves derivative w.r.t. z
def H3(Psi, banded_matrix, mu, nx, ny, nz):
	b = mu*Psi[2:nx, :, :] + (1-2*mu)*Psi[1:nx-1, :, :] + mu*Psi[0:nx-2, :, :]

	Psi_new = np.zeros([nx, ny, nz], dtype = np.complex128)
	Psi_new[1:nx-1, :, :] = solve_banded((1, 1), banded_matrix, b)

	return Psi_new

# Computes the banded matrix
def banded(mu, nx):
	ret_val = np.zeros([3, nx-2], dtype = np.complex128)
	ret_val[0, 1:nx-2] = -mu*np.ones(nx-3)
	ret_val[1, :] = (1 + 2*mu)*np.ones(nx-2)
	ret_val[2, 0:nx-3] = -mu*np.ones(nx-3)

	return ret_val

# To make 2D colormaps of the solution
def plot2D(x, z, prob, t, i):
	plt.close()
	fig = plt.figure()
	plt.axis('square')
	ax = plt.gca()
	probmap = ax.pcolormesh(x, z, prob, cmap = cm.jet, shading = 'auto', vmin = 0.0, vmax = 2.5)
	cbar = plt.colorbar(probmap, orientation='vertical')
	cbar.set_label(r'$|\psi^2|$', fontsize = 14, rotation = 0, labelpad = 20)
	ax.set_title(r't = %.6f' %(t))
	ax.set_xlabel(r'$x$')
	ax.set_ylabel(r'$z$')
	ax.set_xticks([-1, 0, 1])
	ax.set_yticks([-1, 0, 1])
	ax.set_xlim([x[0, 0], x[-1, -1]])
	ax.set_ylim([z[0, 0], z[-1, -1]])
	filename = '/Users/Kishan/Downloads/figures/fig_%03d.png' %(i)
	plt.savefig(filename)
	
	return

# Computes the evolution of Psi with respect to t.
def evolve(nx, ny, nz, m, x, y, z, t, N, skip, arr_len, nprocs, rank):
	dx = x[1] - x[0]
	dz = z[1] - z[0]
	dt = t[1] - t[0]

	x, y, z = np.meshgrid(x, y, z, indexing = 'ij')

	# Parameter mu required to compute banded matrix
	mu = 1j*dt/2/dx**2
	banded_matrix = banded(mu, nx)

	Psi = np.zeros([nx, ny, nz], dtype = np.complex128)
	Psi = init_cond(x, y, z)

	# Computational time and Communication time
	comp_time = 0.0
	comm_time = 0.0

	# 2D
	if ny == 1:
		#plot2D(x[:, 0, :], z[:, 0, :], np.abs(Psi[:, 0, :])**2, t[0], 0)
		#integral[0] = trapz(trapz(np.abs(Psi[:, 0, :])**2, dx = dz), dx = dx)

		# Dividing matrices along z-dimension among different processors
		xv = x[:, :, rank*arr_len:(rank+1)*arr_len]
		yv = y[:, :, rank*arr_len:(rank+1)*arr_len]
		zv = z[:, :, rank*arr_len:(rank+1)*arr_len]
		Psi = Psi[:, :, rank*arr_len:(rank+1)*arr_len]
		nz = arr_len

		#Psi_temp = np.zeros_like(Psi, dtype = np.complex128)

		start = time.time()
		for k in range(1, m):
			start_comp = time.time()
			# Evolution for H0 and H1
			Psi = H0_2D(Psi, xv, yv, zv, nx, ny, nz, dt, N)
			Psi = H1(Psi, banded_matrix, mu, nx, ny, nz)
			end_comp = time.time()
			comp_time = comp_time + end_comp - start_comp

			# Before evolving H2, we need to re-distribute the array along x-axis
			start_comm = time.time()
			for i in range(nprocs):
				Psi[i*arr_len:(i+1)*arr_len, :, :] = np.transpose(Psi[i*arr_len:(i+1)*arr_len, :, :], axes = [2, 1, 0])
			comm.Alltoall(MPI.IN_PLACE, Psi)
			end_comm = time.time()
			comm_time = comm_time + end_comm - start_comm

			start_comp = time.time()
			# Evolution for H3
			Psi = H3(Psi, banded_matrix, mu, nx, ny, nz)
			end_comp = time.time()
			comp_time = comp_time + end_comp - start_comp

			# Again redistributing the array along z-axis for next iteration
			start_comm = time.time()
			for i in range(nprocs):
				Psi[i*arr_len:(i+1)*arr_len, :, :] = np.transpose(Psi[i*arr_len:(i+1)*arr_len, :, :], axes = [2, 1, 0])
			comm.Alltoall(MPI.IN_PLACE, Psi)
			end_comm = time.time()
			comm_time = comm_time + end_comm - start_comm
			'''
			# Plotting and computing integral at specific time steps
			if k % skip == 0:
				if rank == 0:
					Psi_tot = Psi
					for i in range(1, nprocs):
						comm.Recv([Psi_temp, MPI.COMPLEX], source = i, tag = i)
						Psi_tot = np.concatenate((Psi_tot, Psi_temp), axis = 2)

					plot2D(x[:, 0, :], z[:, 0, :], np.abs(Psi_tot[:, 0, :])**2, t[k], int(k/skip))
					integral[int(k/skip)] = trapz(trapz(np.abs(Psi_tot[:, 0, :])**2, dx = dz), dx = dx)

				else:
					comm.Send([Psi, MPI.COMPLEX], dest = 0, tag = rank)
			'''
		stop = time.time()

	# 3D
	else:
		dy = y[1] - y[0]

		# Plotting at a particular value of y-coordinate. # Since we are plotting 2D colormaps, we can't pass a 3D Psi
		y_slice = 16

		#plot2D(x[:, y_slice, :], z[:, y_slice, :], np.abs(Psi[:, y_slice, :])**2, t[0], 0)
		#integral[0] = trapz(trapz(trapz(np.abs(Psi)**2, dx = dz), dx = dy), dx = dx)

		# Dividing matrices along z-dimension among different processors
		xv = x[:, :, rank*arr_len:(rank+1)*arr_len]
		yv = y[:, :, rank*arr_len:(rank+1)*arr_len]
		zv = z[:, :, rank*arr_len:(rank+1)*arr_len]
		Psi = Psi[:, :, rank*arr_len:(rank+1)*arr_len]
		nz = arr_len

		#Psi_temp = np.zeros_like(Psi, dtype = np.complex128)

		start = time.time()
		for k in range(1, m):
			start_comp = time.time()
			# Evolution for H0, H1 and H2
			Psi = H0_3D(Psi, xv, yv, zv, nx, ny, nz, dt, N)
			Psi = H1(Psi, banded_matrix, mu, nx, ny, nz)
			Psi = H2(Psi, banded_matrix, mu, nx, ny, nz)
			end_comp = time.time()
			comp_time = comp_time + end_comp - start_comp

			# Before evolving H2, we need to re-distribute the array along x-axis
			start_comm = time.time()
			for i in range(nprocs):
				Psi[i*arr_len:(i+1)*arr_len, :, :] = np.transpose(Psi[i*arr_len:(i+1)*arr_len, :, :], axes = [2, 1, 0])
			Psi = np.ascontiguousarray(Psi)
			comm.Alltoall(MPI.IN_PLACE, Psi)
			end_comm = time.time()
			comm_time = comm_time + end_comm - start_comm

			start_comp = time.time()
			# Evolution for H3
			Psi = H3(Psi, banded_matrix, mu, nx, ny, nz)
			end_comp = time.time()
			comp_time = comp_time + end_comp - start_comp

			# Again redistributing the array along z-axis for next iteration
			start_comm = time.time()
			for i in range(nprocs):
				Psi[i*arr_len:(i+1)*arr_len, :, :] = np.transpose(Psi[i*arr_len:(i+1)*arr_len, :, :], axes = [2, 1, 0])
			comm.Alltoall(MPI.IN_PLACE, Psi)
			end_comm = time.time()
			comm_time = comm_time + end_comm - start_comm
			'''
			# Plotting and computing integral at specific time steps
			if k % skip == 0:
				if rank == 0:
					Psi_tot = Psi
					for i in range(1, nprocs):
						comm.Recv([Psi_temp, MPI.COMPLEX], source = i, tag = i)
						Psi_tot = np.concatenate((Psi_tot, Psi_temp), axis = 2)

					plot2D(x[:, y_slice, :], z[:, y_slice, :], np.abs(Psi_tot[:, y_slice, :])**2, t[k], int(k/skip))
					integral[int(k/skip)] = trapz(trapz(trapz(np.abs(Psi_tot)**2, dx = dz), dx = dy), dx = dx)

				else:
					comm.Send([Psi, MPI.COMPLEX], dest = 0, tag = rank)
			'''
		stop = time.time()


	if rank == 0:
		print('Time taken:', stop - start)
		print('Comm time:', comm_time)
		print('Comp tiem:', comp_time)

	return


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

arr_len = 64      # Array length per node
nx = arr_len*nprocs
ny = 1
nz = arr_len*nprocs
m = 10001

# Computes the interval at which to plot colormaps.
framerate = 24
vidlen = 6
no_of_frames = framerate*vidlen
skip = int(np.floor((m - 1)/no_of_frames))

# x, y, z, t arrays
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
z = np.linspace(-1, 1, nz)
t = np.linspace(0, 1, m)

# Array to compute integral of |\psi|^2
integral = np.zeros(int(np.floor(m/skip) + 1))
indices = skip*np.linspace(0, len(integral) - 1, len(integral), dtype = np.int)
t_arr = t[indices]

# Parameter to be multiplied with the non-linear term
N = 10.0

evolve(nx, ny, nz, m, x, y, z, t, N, skip, arr_len, nprocs, rank)
'''
# Plotting integral vs t
if rank == 0:
	plt.close()
	plt.plot(t_arr, integral)
	plt.xlabel(r'$t$')
	plt.ylabel(r'$\int|\psi|^2\ \mathregular{d}^2x$')
	plt.savefig('figures/int_par.png', dpi = 300)
'''
