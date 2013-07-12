FC = ifort
FFLAGS = -warn all -check all -check noarg_temp_created
LDFLAGS = -warn all -check all -check noarg_temp_created
FFLAGS_OMP = $(FFLAGS) -openmp
LDFLAGS_OMP = $(LDFLAGS) -openmp
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

module_physics.o : data_precisions.mod data_constants.mod data_cell.mod
class_observables.o : data_precisions.mod data_cell.mod data_monteCarlo.mod
class_neighbourCells.o : data_precisions.mod data_cell.mod data_neighbours.mod module_physics.o
class_mixingPotential.o : data_precisions.mod data_cell.mod data_particles.mod data_potential.mod \
                          data_neighbours.mod module_physics.o class_neighbourCells.o
class_spheres.o : data_precisions.mod data_constants.mod  data_cell.mod module_physics.o \
                  class_neighbourCells.o
class_hardSpheres.o : data_precisions.mod data_cell.mod data_particles.mod data_monteCarlo.mod \
                      data_potential.mod data_neighbours.mod data_distrib.mod module_physics.o \
                      class_observables.o class_neighbourCells.o class_mixingPotential.o \
                      class_spheres.o
class_dipolarSpheres.o : data_precisions.mod data_constants.mod data_cell.mod data_particles.mod \
                         data_monteCarlo.mod data_potential.mod data_neighbours.mod data_distrib.mod \
                         module_physics.o class_observables.o class_neighbourCells.o \
                         class_mixingPotential.o class_spheres.o
module_tools.o : data_precisions.mod data_cell.mod data_particles.mod data_monteCarlo.mod \
                 module_physics.o class_spheres.o

# -----------------------------------------------

module_physics.o : module_physics.f90
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
module_tools.o : module_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@

# -----------------------------------------------

program_canonical.o : program_canonical.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC1) : data.o module_physics.o class_neighbourCells.o class_mixingPotential.o class_spheres.o \
	class_hardSpheres.o class_dipolarSpheres.o class_observables.o class_units.o \
	module_tools.o program_canonical.o
	$(FC) $^ $(LDFLAGS) $(CLDFLAGS) $(LDLIBS) -o $@

# ------------------------------------------------------------------------

program_bunching.o : program_bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC2) : data.o program_bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ------------------------------------------------------------------------

program_distrib.o : program_distrib.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
$(EXEC3) : data.o module_physics.o class_neighbourCells.o class_mixingPotential.o class_spheres.o \
           class_interactingSpheres.o module_tools.o program_distrib.o
	$(FC) $^ $(LDFLAGS_OMP) -o $@

# --------------------------------------------------------------------------------------------------
	
clean :
	rm -f *.o *.mod

mrproper : clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90

remove :
	rm -f *.out *.txt *.shots *.tmp
