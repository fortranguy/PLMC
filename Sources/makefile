FC = ifort
FFLAGS = # -warn all -check all -check noarg_temp_created
LDFLAGS = # -warn all -check all -check noarg_temp_created
#FC = gfortran
#FFLAGS = -W -Wall -fbounds-check
#LDFLAGS = -W -Wall -fbounds-check
EXEC1 = xMC_Canonique
EXEC2 = xBunching
EXEC3 = xDistrib

all: $(EXEC1) $(EXEC2) $(EXEC3)

#----------------------------------------------------------

data_copy.f90: data.f90
	cp data.f90 data_copy.f90

data_copy.o data_constants.mod data_cell.mod data_particles.mod data_mc.mod \
	data_potentiel.mod data_distrib.mod: data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod
	
data_neighbours.o data_neighbours.mod : data_neighbours.f90 data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

# ---------------------------

mod_neighbours.o : data_constants.mod data_cell.mod data_potentiel.mod \
	data_neighbours.mod
mod_tools.o: data_constants.mod data_particles.mod data_mc.mod \
	data_potentiel.mod data_neighbours.mod
mod_physique.o: data_constants.mod data_cell.mod data_particles.mod \
	data_mc.mod data_potentiel.mod data_neighbours.mod mod_neighbours.mod
mc_canonique.o: data_constants.mod data_mc.mod data_distrib.mod\
	mod_physique.mod mod_tools.mod

mod_neighbours.o: mod_neighbours.f90
	$(FC) $(LDFLAGS) -c $<
mod_tools.o: mod_tools.f90
	$(FC) $(LDFLAGS) -c $<
mod_physique.o: mod_physique.f90
	$(FC) $(LDFLAGS) -c $<
mc_canonique.o: mc_canonique.f90
	$(FC) $(LDFLAGS) -c $<

xMC_Canonique: data_copy.o data_neighbours.o mod_neighbours.o mod_physique.o \
	mod_tools.o mc_canonique.o
	$(FC) $(LDFLAGS) $^ -o $@
	
# ---------------------------

bunching.o: data_constants.mod data_mc.mod

bunching.o: bunching.f90
	$(FC) $(FFLAGS) -c $<
	
xBunching: data_copy.o bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ---------------------------

distrib.o: data_constants.mod data_distrib.mod data_particles.mod \
	mod_physique.mod

distrib.o: distrib.f90
	$(FC) $(FFLAGS) -c $<

xDistrib: data_copy.o data_neighbours.o mod_neighbours.o mod_physique.o distrib.o
	$(FC) $^ $(LDFLAGS) -o $@

# ---------------------------------------------------------
	
clean:
	rm -f *.o data_constants.mod data_cell.mod data_particles.mod data_mc.mod \
	data_potentiel.mod data_neighbours.mod
	
mrproper: clean
	rm -f data_copy.f90 *.mod $(EXEC1) $(EXEC2) $(EXEC3) 
