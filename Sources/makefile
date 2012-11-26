FC = ifort # gfortran #
FFLAGS = -warn all -check noarg_temp_created -g # -W -Wall -fbounds-check #
LDFLAGS = -warn all -check noarg_temp_created -g # -W -Wall -fbounds-check #
EXEC1 = xMC_Canonique
EXEC2 = xBunching

all: $(EXEC1) $(EXEC2)

#----------------------------------------------------------

data_copy.f90:
	cp data.f90 data_copy.f90

data_copy.o data_constants.mod data_cell.mod data_neighbours.mod \
    data_particles.mod data_mc.mod data_potentiel.mod: data_copy.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

# ---------------------------

mod_tools.o: data_particles.mod data_mc.mod data_potentiel.mod \
	data_constants.mod
mod_neighbours.o : data_cell.mod data_neighbours.mod
mod_physique.o: data_cell.mod data_particles.mod data_mc.mod \
	data_potentiel.mod data_constants.mod 
mc_canonique.o: data_cell.mod data_particles.mod data_mc.mod \
	data_constants.mod

mod_tools.o: mod_tools.f90
	$(FC) $(LDFLAGS) -c $<
mod_neighbours.o: mod_neighbours.f90
	$(FC) $(LDFLAGS) -c $<
mod_physique.o: mod_physique.f90
	$(FC) $(LDFLAGS) -c $<
mc_canonique.o: mc_canonique.f90
	$(FC) $(LDFLAGS) -c $<

xMC_Canonique: data_copy.o mod_tools.o mod_neighbours.o mod_physique.o \
    mc_canonique.o
	$(FC) $(LDFLAGS) $^ -o $@
	
# ---------------------------

bunching.o: data_mc.mod data_constants.mod

bunching.o: bunching.f90
	$(FC) $(FFLAGS) -c $<
	
xBunching: data_copy.o bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ---------------------------------------------------------
	
clean:
	rm -f *.o data_cell.mod data_neighbours.mod data_particles.mod data_mc.mod \
		data_potentiel.mod data_constants.mod
	
mrproper: clean
	rm -f data_copy.f90 *.mod $(EXEC1) $(EXEC2) 
