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

data.o: data_copy.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_neighbours.o: class_neighbours.f90
	$(FC) -c $(FFLAGS) $< -o $@
class_component.o: class_component.f90
	$(FC) -c $(FFLAGS) $< -o $@
mod_tools.o: mod_tools.f90
	$(FC) -c $(FFLAGS) $< -o $@

# ---------------------------

mc_canonique.o: mc_canonique.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC1) : data.o class_neighbours.o class_component.o mod_tools.o \
	mc_canonique.o
	$(FC) $^ $(LDFLAGS) -o $@
	
# ---------------------------

bunching.o: bunching.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(EXEC2) : data.o bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ---------------------------

distrib.o: distrib.f90
	$(FC) -c $(FFLAGS_OMP) $< -o $@
$(EXEC3) : data.o class_neighbours.o class_component.o mod_tools.o \
	distrib.o
	$(FC) $^ $(LDFLAGS_OMP) -o $@

# ---------------------------------------------------------
	
clean:
	rm -f *.o *.mod
	
mrproper: clean
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) data_copy.f90
