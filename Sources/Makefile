FC = ifort
FFLAGS = -fast -openmp # -warn all -check all -check noarg_temp_created
LDFLAGS = -fast -openmp # -warn all -check all -check noarg_temp_created
FFLAGS_OMP = $(FFLAGS)
LDFLAGS_OMP = $(LDFLAGS)
#FC = gfortran
#FFLAGS = -W -Wall -fbounds-check
#LDFLAGS = -W -Wall -fbounds-check

EXEC1 = xMC_Canonical
EXEC2 = xMC_Bunching
EXEC3 = xMC_Distrib

all: $(EXEC1) $(EXEC2) $(EXEC3)

#---------------------------------------------------------------------------------------------------

data_copy.f90 : data.f90
	cp data.f90 data_copy.f90

data_copy.o data_precisions.mod data_constants.mod data_cell.mod data_particles.mod \
data_monteCarlo.mod data_potential.mod data_neighbours.mod data_monteCarlo_distrib.mod : data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

data.o: data_copy.f90
	$(FC) -c $(FFLAGS) $< -o $@

# ------------------------------------------------------------------------

mod_physics.o : data_precisions.mod data_constants.mod data_cell.mod
class_observables.o : data_precisions.mod data_cell.mod data_monteCarlo.mod
class_neighbourCells.o : data_precisions.mod data_cell.mod data_neighbours.mod mod_physics.o
class_mixingPotential.o : data_precisions.mod data_cell.mod data_particles.mod data_potential.mod \
                          data_neighbours.mod mod_physics.o class_neighbourCells.o
class_spheres.o : data_precisions.mod data_constants.mod  data_cell.mod mod_physics.o \
                  class_neighbourCells.o
class_hardSpheres.o : data_precisions.mod data_cell.mod data_particles.mod data_monteCarlo.mod \
                      data_potential.mod data_neighbours.mod data_distrib.mod mod_physics.o \
                      class_observables.o class_neighbourCells.o class_mixingPotential.o \
                      class_spheres.o
class_dipolarSpheres.o : data_precisions.mod data_constants.mod data_cell.mod data_particles.mod \
                         data_monteCarlo.mod data_potential.mod data_neighbours.mod data_distrib.mod \
                         mod_physics.o class_observables.o class_neighbourCells.o \
                         class_mixingPotential.o class_spheres.o
mod_tools.o : data_precisions.mod data_cell.mod data_particles.mod data_monteCarlo.mod mod_physics.o \
              class_spheres.o

# -----------------------------------------------

mod_physics.o : mod_physics.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
class_observables.o : class_observables.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_neighbourCells.o : class_neighbourCells.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_mixingPotential.o : class_mixingPotential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_spheres.o : class_spheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hardSpheres.o : class_hardSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_interactingSpheres.o : class_interactingSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_dipolarSpheres.o : class_dipolarSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_units.o : class_units.f90
	$(FC) -c $(FFLAGS) $< -o $@
mod_tools.o : mod_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@

# -----------------------------------------------

mc_canonical.o : mc_canonical.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC1) : data.o mod_physics.o class_neighbourCells.o class_mixingPotential.o class_spheres.o \
	class_hardSpheres.o class_dipolarSpheres.o class_observables.o class_units.o \
	mod_tools.o mc_canonical.o
	$(FC) $^ $(LDFLAGS) $(CLDFLAGS) $(LDLIBS) -o $@

# ------------------------------------------------------------------------

mc_bunching.o : mc_bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC2) : data.o mc_bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ------------------------------------------------------------------------

mc_distrib.o : mc_distrib.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
$(EXEC3) : data.o mod_physics.o class_neighbourCells.o class_mixingPotential.o class_spheres.o \
           class_interactingSpheres.o mod_tools.o mc_distrib.o
	$(FC) $^ $(LDFLAGS_OMP) -o $@

# --------------------------------------------------------------------------------------------------
	
clean :
	rm -f *.o *.mod

mrproper : clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90

remove :
	rm -f *.out *.txt *.shots *.tmp
