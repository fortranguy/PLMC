# Monte-Carlo simulation Makefile

ifeq ($(compiler), gfortran)
	FC = gfortran
	# default: debug
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
	# default: debug
	FFLAGS = -warn all -check all -check noarg_temp_created -fltconsistency -g -traceback
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
EXEC3 = mc_histogram
EXEC4 = mc_distribution

all: $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4)
	# Report
	git log | head -n3 > make_report.txt
	git diff --color-words >> make_report.txt
	echo "FC = ${FC}" >> make_report.txt
	echo "FFLAGS = ${FFLAGS}" >> make_report.txt
	echo "LDFLAGS = ${LDFLAGS}" >> make_report.txt

#---------------------------------------------------------------------------------------------------

data.o: data.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_types_micro.o: module_types_micro.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_physics_micro.o: module_physics_micro.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_small_move.o: class_small_move.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_small_rotation.o: class_small_rotation.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hard_spheres.o: class_hard_spheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_dipolar_spheres.o: class_dipolar_spheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_neighbour_cells.o: class_neighbour_cells.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hard_spheres_potential.o: class_hard_spheres_potential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_mixing_potential.o: class_mixing_potential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_observables.o: class_observables.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_units.o: class_units.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_physics_macro.o: module_physics_macro.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_algorithms.o: module_algorithms.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_write.o: module_write.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_monte_carlo_arguments.o: module_monte_carlo_arguments.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_physical_system.o: class_physical_system.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_monte_carlo.o: program_monte_carlo.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_bunching.o: program_bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_histogram.o: program_histogram.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_distribution.o: program_distribution.f90
	$(FC) -c $(FFLAGS) $< -o $@
	
# ------------------------------------------------------------------------

$(EXEC1): data.o \
          module_types_micro.o module_physics_micro.o \
          class_small_move.o class_small_rotation.o\
          class_hard_spheres.o class_dipolar_spheres.o \
          class_neighbour_cells.o \
          class_hard_spheres_potential.o \
          class_mixing_potential.o \
          class_observables.o class_units.o \
          module_physics_macro.o module_algorithms.o module_write.o module_monte_carlo_arguments.o\
          class_physical_system.o \
          program_monte_carlo.o
	$(FC) $^ $(LDFLAGS) -o $@
	
$(EXEC2): data.o program_bunching.o
	$(FC) $^ $(LDFLAGS) -o $@
	
$(EXEC3): data.o program_histogram.o
	$(FC) $^ $(LDFLAGS) -o $@
	
$(EXEC4): data.o module_physics_micro.o program_distribution.o
	$(FC) $^ $(LDFLAGS) -o $@

# --------------------------------------------------------------------------------------------------
.PHONY: clean mrproper remove

clean:
	rm -f *.o *.mod

mrproper: clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4)

remove:
	rm -f *.out *.txt *.shots *.tmp
