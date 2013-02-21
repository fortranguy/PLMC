FC = ifort
FFLAGS = #-warn all -check all -check noarg_temp_created
LDFLAGS = #-warn all -check all -check noarg_temp_created
FFLAGS_OMP = $(FFLAGS) -openmp 
LDFLAGS_OMP = $(LDFLAGS) -openmp
#FC = gfortran
#FFLAGS = -W -Wall -fbounds-check
#LDFLAGS = -W -Wall -fbounds-check
EXEC1 = xMC_Canonique
EXEC2 = xBunching
EXEC3 = xDistrib

all: $(EXEC1) $(EXEC2) $(EXEC3)

#----------------------------------------------------------

data_copy.f90 : data.f90	
	cp data.f90 data_copy.f90

# ---------------------------

$(EXEC1) : data.f90 data_neighbours.f90 mod_neighbours.f90 mod_physique.f90 \
	mod_tools.f90 class_particle.f90 obj_particles.f90  mc_canonique.f90
	$(FC) $(LDFLAGS) $^ -o $@
	
# ---------------------------

$(EXEC2) : data_copy.f90 bunching.f90
	$(FC) $^ $(LDFLAGS) -o $@

# ---------------------------

$(EXEC3) : data.f90 data_neighbours.f90 mod_neighbours.f90 mod_physique.f90 \
	distrib.f90
	$(FC) $^ $(LDFLAGS_OMP) -o $@

# ---------------------------------------------------------
	
clean:
	rm -f *.o *.mod
	
mrproper: clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90
