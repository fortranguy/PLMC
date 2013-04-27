FC = ifort
FFLAGS = #-warn all -check all -check noarg_temp_created
LDFLAGS = #-warn all -check all -check noarg_temp_created
FFLAGS_OMP = $(FFLAGS) -openmp 
LDFLAGS_OMP = $(LDFLAGS) -openmp
#FC = gfortran
#FFLAGS = -W -Wall -fbounds-check
#LDFLAGS = -W -Wall -fbounds-check
EXEC1 = xMC_Canonical
EXEC2 = xBunching
EXEC3 = xDistrib

all: $(EXEC1) $(EXEC2) $(EXEC3)

#---------------------------------------------------------------------------------------------------

data_copy.f90 : data.f90
	cp data.f90 data_copy.f90
	
data_copy.o data_constants.mod data_cell.mod data_particles.mod data_mc.mod data_potentiel.mod \
	data_neighbours.mod data_distrib.mod : data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

data.o: data_copy.f90
	$(FC) -c $(FFLAGS) $< -o $@
	
# ------------------------------------------------------------------------

mod_physics.o : data_constants.mod data_cell.mod
class_neighbours.o : data_constants.mod data_cell.mod data_neighbours.mod
class_mixingPotential.o : data_constants.mod data_particles.mod data_mc.mod data_potentiel.mod \
	data_neighbours.mod \
	mod_physics.o class_neighbours.o
class_spheres.o : data_constants.mod data_cell.mod data_mc.mod data_neighbours.mod \
	mod_physics.o class_neighbours.o
class_hardSpheres.o : data_constants.mod data_cell.mod data_particles.mod data_potentiel.mod \
	data_mc.mod data_neighbours.mod \
	mod_physics.o class_neighbours.o class_mixingPotential.o class_spheres.o
class_interactingSpheres.o : data_constants.mod data_cell.mod \
	data_particles.mod data_potentiel.mod data_mc.mod data_neighbours.mod \
	mod_physics.o class_neighbours.o class_mixingPotential.o class_spheres.o
class_observables.o : data_constants.mod data_mc.mod
mod_tools.o : data_constants.mod data_particles.mod data_mc.mod \
	mod_physics.o class_spheres.o

# -----------------------------------------------

mod_physics.o : mod_physics.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
class_neighbours.o : class_neighbours.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_mixingPotential.o : class_mixingPotential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_spheres.o : class_spheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hardSpheres.o : class_hardSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_interactingSpheres.o : class_interactingSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_observables.o : class_observables.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_units.o : class_units.f90
	$(FC) -c $(FFLAGS) $< -o $@
mod_tools.o : mod_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@

# -----------------------------------------------

mc_canonical.o : mc_canonical.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC1) : data.o mod_physics.o class_neighbours.o class_mixingPotential.o class_spheres.o \
	class_hardSpheres.o class_interactingSpheres.o  class_observables.o class_units.o \
	mod_tools.o mc_canonical.o
	$(FC) $^ $(LDFLAGS) -o $@
	
# ------------------------------------------------------------------------

bunching.o : bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC2) : data.o bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ------------------------------------------------------------------------

distrib.o : distrib.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
$(EXEC3) : data.o mod_physics.o class_neighbours.o class_spheres.o \
	class_interactingSpheres.o mod_tools.o distrib.o
	$(FC) $^ $(LDFLAGS_OMP) -o $@

# --------------------------------------------------------------------------------------------------
	
clean :
	rm -f *.o *.mod
	
mrproper : clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90
