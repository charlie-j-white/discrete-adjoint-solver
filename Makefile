#**********************************************************************#
#                          Charlie Anderson                            #
#                  University of Bristol, March 2019                   #
#**********************************************************************#


# define fortran compiler options
FC =gfortran
FLIBS =-llapack
FFLAGS =-O3 -Wall -fbounds-check


# define c compiler options
CC =gcc
CFLAGS =-lm -O3 -Wall


# specify directory for source code files
SDIR =SRC


# name which object files to use
_OBJS = wrapper.o\
       main.o\
       meshing.o\
       split_fwd.o\
       split_rev.o\
       debug_cell.o\
       postprocess.o\
       initialise.o\
       resid.o\
       boundaries.o\
       update.o\
       pressure.o\
       fg_vector.o\
       itinfo.o\
       jst_calcs.o\
       timestep.o\
       cost_is.o\
       aresid.o\
       aMcost_is_d.o\
       aMresid_d.o\
       adjoint.o\
       colour.o\
       lodepng.o


# specify directory for object files and put them there
ODIR =OBJ
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


# first rule, compile all *.f files
$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) -c -o $@ $< $(FFLAGS)


# second rule, compile all *.c files
$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS)


# main rule, performs linking
main: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS) $(FLIBS)





# clean rules, `refreshes' the build by deleting the necessary files
.PHONY: clean
clean:
	rm $(OBJS) main *.plt *.dat *.png PROCESS


.PHONY: cleanup
cleanup:
	rm $(OBJS) main











