FC = ifort
FFLAGS = -warn all -check all -check noarg_temp_created # -fast -openmp #
LDFLAGS = -warn all -check all -check noarg_temp_created # -fast -openmp #
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

data_copy.o data_precisions.mod data_constants.mod data_box.mod data_particles.mod \
data_monteCarlo.mod data_potential.mod data_neighbourCells.mod data_distribution.mod : data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

data.o: data_copy.f90
	$(FC) -c $(FFLAGS) $< -o $@

# ------------------------------------------------------------------------

module_physics.o : data_precisions.mod data_constants.mod data_box.mod
class_observables.o : data_precisions.mod data_box.mod data_monteCarlo.mod
class_neighbourCells.o : data_precisions.mod data_box.mod data_neighbourCells.mod \
                         module_physics.o
class_spheres.o : data_precisions.mod data_constants.mod data_box.mod \
                  module_physics.o \
                  class_neighbourCells.o
class_mixingPotential.o : data_precisions.mod data_box.mod data_particles.mod data_potential.mod \
                          data_neighbourCells.mod \
                          module_physics.o \
                          class_neighbourCells.o class_spheres.o
class_dipolarSpheres.o : data_precisions.mod data_constants.mod data_box.mod data_particles.mod \
                         data_monteCarlo.mod data_potential.mod data_neighbourCells.mod \
                         data_distribution.mod \
                         module_physics.o \
                         class_observables.o class_neighbourCells.o class_mixingPotential.o \
                         class_spheres.o
class_hardSpheres.o : data_precisions.mod data_box.mod data_particles.mod data_monteCarlo.mod \
                      data_neighbourCells.mod data_distribution.mod \
                      module_physics.o \
                      class_observables.o class_neighbourCells.o class_mixingPotential.o \
                      class_spheres.o
classe_spheres_change.o : class_spheres.o class_dipolarSpheres.o class_hardSpheres.o
module_tools.o : data_precisions.mod data_box.mod data_monteCarlo.mod \
                 module_physics.o \
                 class_spheres.o class_dipolarSpheres.o

# -----------------------------------------------

module_physics.o : module_physics.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
class_observables.o : class_observables.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_neighbourCells.o : class_neighbourCells.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_spheres.o : class_spheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_mixingPotential.o : class_mixingPotential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hardSpheres.o : class_hardSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_interactingSpheres.o : class_interactingSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_dipolarSpheres.o : class_dipolarSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
classe_spheres_change.o : classe_spheres_change.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_units.o : class_units.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_tools.o : module_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@

# -----------------------------------------------

program_monte-carlo.o : program_monte-carlo.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC1) : data.o module_physics.o class_neighbourCells.o class_spheres.o class_mixingPotential.o \
	class_hardSpheres.o class_dipolarSpheres.o class_interactingSpheres.o class_observables.o \
	class_units.o module_tools.o program_monte-carlo.o
	$(FC) $^ $(LDFLAGS) $(CLDFLAGS) $(LDLIBS) -o $@

# ------------------------------------------------------------------------

program_bunching.o : program_bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC2) : data.o program_bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ------------------------------------------------------------------------

program_distribution.o : program_distribution.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
$(EXEC3) : data.o module_physics.o class_neighbourCells.o class_spheres.o class_mixingPotential.o \
           class_hardSpheres.o module_tools.o program_distribution.o
	$(FC) $^ $(LDFLAGS_OMP) -o $@

# --------------------------------------------------------------------------------------------------
	
clean :
	rm -f *.o *.mod

mrproper : clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90

remove :
	rm -f *.out *.txt *.shots *.tmp
