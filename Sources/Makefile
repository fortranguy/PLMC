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
	git diff >> make_report.txt
	echo "FC = ${FC}" >> make_report.txt
	echo "FFLAGS = ${FFLAGS}" >> make_report.txt
	echo "LDFLAGS = ${LDFLAGS}" >> make_report.txt

#---------------------------------------------------------------------------------------------------

data.o: data.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_types.o: module_types.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_physics_micro.o: module_physics_micro.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_neighbourCells.o: class_neighbourCells.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_hardSpheres.o: class_hardSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_dipolarSpheres.o: class_dipolarSpheres.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_mixingPotential.o: class_mixingPotential.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_observables.o: class_observables.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_units.o: class_units.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_physics_macro.o: module_physics_macro.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_algorithms.o: module_algorithms.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_tools.o: module_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@
module_distribution.o: module_distribution.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_monte-carlo.o: program_monte-carlo.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_bunching.o: program_bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_histogram.o: program_histogram.f90
	$(FC) -c $(FFLAGS) $< -o $@
program_distribution.o: program_distribution.f90
	$(FC) -c $(FFLAGS) $< -o $@
	
# ------------------------------------------------------------------------

$(EXEC1): data.o \
          module_types.o module_physics_micro.o \
          class_neighbourCells.o \
          class_hardSpheres.o class_dipolarSpheres.o \
          class_mixingPotential.o \
          class_observables.o class_units.o \
          module_physics_macro.o module_algorithms.o module_tools.o \
          program_monte-carlo.o
	$(FC) $^ $(LDFLAGS) -o $@
	
$(EXEC2): data.o program_bunching.o
	$(FC) $^ $(LDFLAGS) -o $@
	
$(EXEC3): data.o program_histogram.o
	$(FC) $^ $(LDFLAGS) -o $@
	
$(EXEC4): data.o module_physics_micro.o module_distribution.o program_distribution.o
	$(FC) $^ $(LDFLAGS) -o $@

# --------------------------------------------------------------------------------------------------
.PHONY: clean mrproper remove

clean:
	rm -f *.o *.mod

mrproper: clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4)

remove:
	rm -f *.out *.txt *.shots *.tmp
