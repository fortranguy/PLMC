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

FFLAGS += -I${HOME}/json-fortran/include
LDFLAGS += -L${HOME}/json-fortran/lib -ljson

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

$(EXEC1): data.o \
          module_types_micro.o module_physics_micro.o module_data.o\
          class_small_move.o class_small_rotation.o\
          class_hard_spheres.o \
          class_neighbour_cells.o \
          class_hard_spheres_potential_energy.o class_ewald_summation_real.o class_ewald_summation_reci.o\
          class_ewald_summation_self.o class_ewald_summation_bound.o\
          class_mixing_potential_energy.o \
          class_hard_spheres_observables.o class_hard_spheres_units.o \
          module_types_macro.o \
          module_physics_macro.o module_algorithms.o module_write.o module_monte_carlo_arguments.o\
          class_physical_system.o \
          program_monte_carlo.o
	$(FC) $^ -o $@ $(LDFLAGS)
	
$(EXEC2): data.o program_bunching.o
	$(FC) $^ -o $@ $(LDFLAGS)
	
$(EXEC3): data.o program_histogram.o
	$(FC) $^ -o $@ $(LDFLAGS)
	
$(EXEC4): data.o module_physics_micro.o program_distribution.o
	$(FC) $^ -o $@ $(LDFLAGS)
	
# ------------------------------------------------------------------------

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# --------------------------------------------------------------------------------------------------
.PHONY: clean mrproper remove

clean:
	rm -f *.o *.mod

mrproper: clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4)

remove:
	rm -f *.out *.txt *.shots *.tmp
