# Monte-Carlo simulation Makefile

ifeq ($(compiler), gfortran)

	FC = gfortran
	# default : debug
	FFLAGS = -W -Wall -fbounds-check -g
	ifeq ($(mode), release)
		FFLAGS = -O2 -g
	endif
	ifeq ($(mode), fast)
		FFLAGS = -Ofast
	endif
	ifeq ($(omp), true)
        FFLAGS += -fopenmp
        LDFLAGS = -fopenmp
    endif
	
else

	FC = ifort
	# default : debug
	FFLAGS = -warn all -check all -check noarg_temp_created -fltconsistency -g
	ifeq ($(mode), release)
		FFLAGS =  -fltconsistency -g
	endif
	ifeq ($(mode), fast)
		FFLAGS = -fast
	endif
	ifeq ($(omp), true)
        FFLAGS += -openmp
        LDFLAGS = -openmp
    endif
	
endif

# Default
compiler = ifort
mode = debug

EXEC1 = mc_cano_bulk
EXEC2 = mc_bunching
EXEC3 = mc_distribution

all: $(EXEC1) $(EXEC2) $(EXEC3)
	# Report
	git log | head -n3 > make_report.txt
	git diff >> make_report.txt
	echo "FC = ${FC}" >> make_report.txt
	echo "FFLAGS = ${FFLAGS}" >> make_report.txt
	echo "LDFLAGS = ${LDFLAGS}" >> make_report.txt

#---------------------------------------------------------------------------------------------------

data_copy.f90 : data.f90
	cp data.f90 data_copy.f90
	chmod ugo-w data_copy.f90

data_copy.o data_precisions.mod data_constants.mod data_box.mod data_particles.mod \
data_monteCarlo.mod data_potential.mod data_neighbourCells.mod data_distribution.mod : data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

data.o: data_copy.f90
	$(FC) -c $(FFLAGS) $< -o $@

# ------------------------------------------------------------------------

module_physics.o : data_precisions.mod data_constants.mod data_box.mod
class_neighbourCells.o : data_precisions.mod data_box.mod data_neighbourCells.mod \
                         module_physics.o
class_hardSpheres.o : data_precisions.mod data_constants.mod data_box.mod data_particles.mod \
                      data_monteCarlo.mod data_neighbourCells.mod data_distribution.mod \
                      module_physics.o \
                      class_neighbourCells.o
class_interactingSpheres.o : data_precisions.mod data_box.mod data_particles.mod data_potential.mod \
                             data_monteCarlo.mod data_neighbourCells.mod data_distribution.mod \
                             module_physics.o \
                             class_neighbourCells.o class_hardSpheres.o
class_dipolarSpheres.o : data_precisions.mod data_constants.mod data_box.mod data_particles.mod \
                         data_monteCarlo.mod data_potential.mod data_neighbourCells.mod \
                         data_distribution.mod \
                         module_physics.o \
                         class_neighbourCells.o class_hardSpheres.o
class_mixingPotential.o : data_precisions.mod data_box.mod data_particles.mod data_potential.mod \
                          data_neighbourCells.mod \
                          module_physics.o \
                          class_neighbourCells.o class_hardSpheres.o
class_observables.o : data_precisions.mod data_box.mod data_monteCarlo.mod
module_algorithms.o : data_precisions.mod data_box.mod data_monteCarlo.mod \
                      module_physics.o \
                      class_observables.mod class_hardSpheres.o class_interactingSpheres.o \
                      class_dipolarSpheres.o class_mixingPotential.o
module_tools.o : data_precisions.mod data_box.mod data_monteCarlo.mod \
                 module_physics.o \
                 class_hardSpheres.o class_interactingSpheres.o class_dipolarSpheres.o \
                 class_mixingPotential.o class_observables.o class_units.o
module_distribution.o : data_precisions.mod data_constants.mod data_distribution.mod

# -----------------------------------------------

module_physics.o : module_physics.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_neighbourCells.o : class_neighbourCells.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hardSpheres.o : class_hardSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_interactingSpheres.o : class_interactingSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_dipolarSpheres.o : class_dipolarSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_mixingPotential.o : class_mixingPotential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_observables.o : class_observables.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_units.o : class_units.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_algorithms.o : module_algorithms.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_tools.o : module_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_distribution.o : module_distribution.f90
	$(FC) -c $(FFLAGS) $< -o $@

# -----------------------------------------------

program_monte-carlo.o : program_monte-carlo.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC1) : data.o module_physics.o class_neighbourCells.o class_hardSpheres.o \
           class_interactingSpheres.o class_dipolarSpheres.o class_mixingPotential.o \
           class_observables.o class_units.o \
           module_algorithms.o module_tools.o program_monte-carlo.o
	$(FC) $^ $(LDFLAGS) -o $@

# ------------------------------------------------------------------------

program_bunching.o : program_bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC2) : data.o program_bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ------------------------------------------------------------------------

program_distribution.o : program_distribution.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC3) : data.o module_physics.o module_distribution.o program_distribution.o
	$(FC) $^ $(LDFLAGS) -o $@

# --------------------------------------------------------------------------------------------------
.PHONY : clean mrproper remove

clean :
	rm -f *.o *.mod

mrproper : clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90

remove :
	rm -f *.out *.txt *.shots *.tmp
