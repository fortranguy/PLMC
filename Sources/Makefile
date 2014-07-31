# Monte-Carlo simulation Makefile

ifeq ($(compiler), gfortran)
    FC = gfortran
    CPPC = g++
    # default: debug
    FFLAGS = -W -Wall -fbounds-check -g
    CPPFLAGS = -W -Wall -g
    ifeq ($(mode), release)
        FFLAGS = -O2 -g
        CPPFLAGS = -O2
    endif
    ifeq ($(mode), fast)
        FFLAGS = -Ofast
    endif
    FFLAGS += -I${HOME}/usr/json-fortran_gfortran/include
    LDFLAGS += -L${HOME}/usr/json-fortran_gfortran/lib -ljsonfortran
else
    FC = ifort
    CPPC = icpc
    # default: debug
    FFLAGS = -warn all -check all -check noarg_temp_created -fltconsistency -g -traceback
    CPPFLAGS = -W -Wall -g
    ifeq ($(mode), release)
        FFLAGS =  -fltconsistency -g
        CPPFLAGS = -O2
    endif
    ifeq ($(mode), fast)
        FFLAGS = -fast
        CPPFLAGS = -O3
    endif
    FFLAGS += -I${HOME}/usr/json-fortran_intel/include
    LDFLAGS += -L${HOME}/usr/json-fortran_intel/lib -ljsonfortran
endif

CPPFLAGS += -I${HOME}/usr/lib/boost_1_55_0
LDLIBS = -lstdc++

# Default
compiler = ifort
mode = debug

all: mc_cano_bulk mc_bunching mc_histogram mc_distribution mc_between_distribution mc_clusters
	# Report
	git log | head -n3 > make_report.txt
	git diff --color-words >> make_report.txt
	echo "FC = ${FC}" >> make_report.txt
	echo "FFLAGS = ${FFLAGS}" >> make_report.txt
	echo "LDFLAGS = ${LDFLAGS}" >> make_report.txt

#---------------------------------------------------------------------------------------------------

mc_cano_bulk: data_precisions.o data_constants.o data_box.o data_neighbour_cells.o data_write.o \
              module_data.o module_types_micro.o module_physics_micro.o \
              class_hard_spheres.o \
              class_small_move.o class_small_rotation.o \
              class_neighbour_cells.o class_hard_spheres_potential.o \
              class_ewald_summation_real.o class_ewald_summation_reci.o \
              class_ewald_summation_self.o class_ewald_summation_bound.o \
              module_types_macro.o \
              class_discrete_observable.o class_hard_spheres_observables.o class_hard_spheres_units.o \
              module_physics_macro.o module_algorithms.o module_write.o \
              class_physical_system.o module_arguments.o module_arguments_monte_carlo.o \
              program_monte_carlo.o
	$(FC) $^ -o $@ $(LDFLAGS)

mc_bunching: module_data.o \
             module_arguments.o \
             program_bunching.o
	$(FC) $^ -o $@ $(LDFLAGS)

mc_histogram: module_data.o module_arguments.o \
              program_histogram.o
	$(FC) $^ -o $@ $(LDFLAGS)

mc_distribution: data_constants.o data_box.o \
                 module_data.o module_physics_micro.o \
                 module_arguments.o \
                 program_distribution.o
	$(FC) $^ -o $@ $(LDFLAGS)
	
mc_between_distribution: data_constants.o data_box.o \
                         module_data.o module_physics_micro.o \
                         module_arguments.o \
                         program_between_distribution.o
	$(FC) $^ -o $@ $(LDFLAGS)

mc_clusters: data_precisions.o data_constants.o data_box.o \
             module_data.o module_physics_micro.o connected_components.o module_clusters.o \
             module_arguments.o \
             program_clusters.o
	$(FC) $^ $(LDFLAGS) -o $@ $(LDLIBS)

# ------------------------------------------------------------------------

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.cpp
	$(CPPC) -c $(CPPFLAGS) $< -o $@

# --------------------------------------------------------------------------------------------------
.PHONY: clean mrproper remove

clean:
	rm -f *.o *.mod

mrproper: clean
	rm -f mc_cano_bulk mc_bunching mc_histogram mc_distribution mc_between_distribution mc_clusters

remove:
	rm -f *.out *.txt *.shots *.tmp
