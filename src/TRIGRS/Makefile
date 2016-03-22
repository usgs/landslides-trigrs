TRG = trg
PRG = prg
TPX = tpx

OBJT90 = trigrs.o flux.o prpijz.o svxmdv.o svijz.o dzero_brac.o
OBJT95 = grids.o input_vars.o model_vars.o dsimps.o input_file_defs.o iverson.o pstpi.o satfin.o savage.o steady.o trini.o unsinf.o ivestp.o pstpf.o rnoff.o satinf.o smallt.o svgstp.o unsfin.o unsth.o ssizgrd.o svlist.o
OBJT77 = calerf.o dbsct.o derfc.o irdgrd.o irdswm.o isvgrd.o roots.o srdgrd.o srdswm.o ssvgrd.o

OBJP90 = trigrs_p.o partial_p.o flux.o flux_p.o prpijz.o svxmdv.o svijz.o dzero_brac.o srdgrd_p.o irdgrd_p.o
OBJP95 = modules_p.o grids.o input_vars.o model_vars.o dsimps.o input_file_defs.o iverson.o pstpi.o pstpi_p.o satfin.o satfin_p.o savage.o steady.o trini.o trini_p.o unsinf.o unsinf_p.o ivestp.o ivestp_p.o pstpf.o pstpf_p.o rnoff.o satinf.o satinf_p.o smallt.o svgstp.o svgstp_p.o unsfin.o unsfin_p.o unsth.o unsth_p.o ssizgrd.o ssizgrd_p.o svlist.o rnoff_p.o steady_p.o
OBJP77 = calerf.o dbsct.o derfc.o irdgrd.o irdswm.o irdswm_p.o isvgrd.o roots.o srdgrd.o srdswm.o srdswm_p.o ssvgrd.o

OBJX90 = tpindx.o nxtcel.o
OBJX95 = ssizgrd.o
OBJX77 = isvgrd.o mpfldr.o rdflodir.o sindex.o slofac.o srdgrd1.o

LIBS =

CC = gcc -O3
CCFLAGS = -lgsl -lgslcblas -lm
FC = f95 -w -O3 
FFLAGS =
F90 = f95 -w -O3
MPIF90 = mpif90 -w -O3
F90FLAGS = -w -O3
LDFLAGS = -w -O3

all: $(TRG) $(PRG) $(TPX)
#-----------------------------------------------------------------

$(TRG): $(OBJT95) $(OBJT90) $(OBJT77) 
	$(F90) $(CCLIBS) $(LDFLAGS) -o $@ $(OBJT95) $(OBJT90) $(OBJT77) $(CCFLAGS) $(LIBS)

$(PRG): $(OBJP95) $(OBJP90) $(OBJP77) 
	$(MPIF90) $(CCLIBS) $(LDFLAGS) -o $@ $(OBJP95) $(OBJP90) $(OBJP77) $(CCFLAGS) $(LIBS)

$(TPX): $(OBJX95) $(OBJX90) $(OBJX77)
	$(F90) $(CCLIBS) $(LDFLAGS) -o $@ $(OBJX95) $(OBJX90) $(OBJX77) $(CCFLAGS) $(LIBS)

#-----------------------------------------------------------------
clean:
	rm -f $(TRG) $(TPX) $(PRG)
	rm -rf $(OBJT95) $(OBJT90) $(OBJT77) $(OBJP95) $(OBJP90) $(OBJP77)
	rm -rf $(OBJX95) $(OBJX90) $(OBJX77)
	rm -rf *.mod *.exe *.stackdump

.SUFFIXES: $(SUFFIXES) .f90 .f .c .f95

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

.f.o:
	$(MPIF90) $(F90FLAGS) -c $<

.c.o:
	$(CC) $(CCINCLUDE) -c -w $<

.f95.o:
	$(MPIF90) $(F90FLAGS) -c $<

