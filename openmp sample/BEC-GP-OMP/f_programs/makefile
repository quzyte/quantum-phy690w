#
# BEC-GP-OMP codes are developed and (c)opyright-ed by:
#
# Luis E. Young-S., Sadhan K. Adhikari
# (UNESP - Sao Paulo State University, Brazil)
#
# Paulsamy Muruganandam
# (Bharathidasan University, Tamil Nadu, India)
#
# Dusan Vudragovic, Antun Balaz
# (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
#
# Public use and modification of this code are allowed provided that the
# following three papers are cited:
#
# [1] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209. 
# [2] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
# [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
#
# The authors would be grateful for all information and/or comments
# regarding the use of the code.
#

FC = ifort
FFLAGS = -O3 -parallel -par-report -V -mcmodel medium -shared-intel 

ifeq ($(compiler), gfort)
FC = gfortran
FFLAGS =  -O3 -fopenmp 
endif

all: imag1d real1d imag2d real2d imagcir realcir imagsph realsph imagaxi realaxi imag3d real3d

help: ../readme.txt
	less $^

imag1d:
	$(FC) $(FFLAGS) -c src/imag1d.f90 -o imag1d.o
	$(FC) $(FFLAGS) imag1d.o -o imag1d 
	rm -rf *.o *.mod

real1d:
	$(FC) $(FFLAGS) -c src/real1d.f90 -o real1d.o
	$(FC) $(FFLAGS) real1d.o -o real1d 
	rm -rf *.o *.mod
	
imagcir:
	$(FC) $(FFLAGS) -c src/imagcir.f90 -o imagcir.o
	$(FC) $(FFLAGS) imagcir.o -o imagcir 
	rm -rf *.o *.mod

realcir:
	$(FC) $(FFLAGS) -c src/realcir.f90 -o realcir.o
	$(FC) $(FFLAGS) realcir.o -o realcir 
	rm -rf *.o *.mod	

imagsph:
	$(FC) $(FFLAGS) -c src/imagsph.f90 -o imagsph.o
	$(FC) $(FFLAGS) imagsph.o -o imagsph 
	rm -rf *.o *.mod	

realsph:
	$(FC) $(FFLAGS) -c src/realsph.f90 -o realsph.o
	$(FC) $(FFLAGS) realsph.o -o realsph 
	rm -rf *.o *.mod	

imag2d:
ifeq ($(compiler), gfort)
	$(FC) $(FFLAGS) -fno-automatic -c src/imag2d.f90 -o imag2d.o
	$(FC) $(FFLAGS) -fno-automatic  imag2d.o -o imag2d 
	rm -rf *.o *.mod
else
	$(FC) $(FFLAGS) -c src/imag2d.f90 -o imag2d.o
	$(FC) $(FFLAGS) imag2d.o -o imag2d 
	rm -rf *.o *.mod
endif

real2d:
	$(FC) $(FFLAGS) -c src/real2d.f90 -o real2d.o
	$(FC) $(FFLAGS) real2d.o -o real2d 
	rm -rf *.o *.mod

imagaxi:
	$(FC) $(FFLAGS) -c src/imagaxi.f90 -o imagaxi.o
	$(FC) $(FFLAGS) imagaxi.o -o imagaxi 
	rm -rf *.o *.mod

realaxi:
	$(FC) $(FFLAGS) -c src/realaxi.f90 -o realaxi.o
	$(FC) $(FFLAGS) realaxi.o -o realaxi 
	rm -rf *.o *.mod
	
imag3d:
ifeq ($(compiler), gfort)
	$(FC) $(FFLAGS) -fno-automatic -c src/imag3d.f90 -o imag3d.o
	$(FC) $(FFLAGS) -fno-automatic imag3d.o -o imag3d 
	rm -rf *.o *.mod
else
	$(FC) $(FFLAGS) -c src/imag3d.f90 -o imag3d.o
	$(FC) $(FFLAGS) imag3d.o -o imag3d 
	rm -rf *.o *.mod
endif

real3d:
ifeq ($(compiler), gfort)
	$(FC) $(FFLAGS) -fno-automatic -c src/real3d.f90 -o real3d.o
	$(FC) $(FFLAGS) -fno-automatic real3d.o -o real3d 
	rm -rf *.o *.mod
else
	$(FC) $(FFLAGS) -c src/real3d.f90 -o real3d.o
	$(FC) $(FFLAGS) real3d.o -o real3d 
	rm -rf *.o *.mod
endif

full-clean: clean # Cleans all the output files also
	rm -rf *.txt 

clean:
	rm -rf *~ *.o imag1d real1d imag2d real2d imagcir realcir imagsph realsph imagaxi realaxi imag3d real3d
