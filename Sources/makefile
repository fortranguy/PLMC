FC = ifort # gfortran #
FFLAGS =  # -warn all -check all # -W -Wall -fbounds-check #
LDFLAGS = # -warn all -check all # -W -Wall -fbounds-check #
EXEC1 = xMC_Canonique
EXEC2 = xBunching

all: $(EXEC1) $(EXEC2)

#----------------------------------------------------------

data.o data_constants.mod data_cell.mod data_particles.mod data_mc.mod \
	data_potentiel.mod: data.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

# ---------------------------

mod_tools.o: data_particles.mod data_mc.mod data_potentiel.mod \
	data_constants.mod
mod_physique.o: data_cell.mod data_particles.mod data_mc.mod \
	data_potentiel.mod data_constants.mod 
mc_canonique.o: data_cell.mod data_particles.mod data_mc.mod \
	data_constants.mod

mod_tools.o: mod_tools.f90
	$(FC) $(LDFLAGS) -c $<
mod_physique.o: mod_physique.f90
	$(FC) $(LDFLAGS) -c $<
mc_canonique.o: mc_canonique.f90
	$(FC) $(LDFLAGS) -c $<

xMC_Canonique: data.o mod_tools.o mod_physique.o mc_canonique.o
	$(FC) $(LDFLAGS) $^ -o $@
	
# ---------------------------

bunching.o: data_mc.mod data_constants.mod

bunching.o: bunching.f90
	$(FC) $(FFLAGS) -c $<
	
xBunching: data.o bunching.o
	$(FC) $^ $(LDFLAGS) -o $@

# ---------------------------------------------------------
	
clean:
	rm -rf *.o data_cell.mod data_particles.mod data_mc.mod \
		data_potentiel.mod data_constants.mod
	
mrproper: clean
	rm -rf *.mod
