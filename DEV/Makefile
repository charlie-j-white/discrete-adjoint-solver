#**********************************************************************#
#                          Charlie Anderson                            #
#                  University of Bristol, March 2019                   #
#**********************************************************************#


# define fortran compiler options
FC =gfortran
FLIBS =
FFLAGS =-O3 -Wall -fbounds-check




# name which object files to use
OBJS = memory.o


# compile all *.f files
%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)



# main rule, performs linking
memory: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS) $(FLIBS)





# clean rules, `refreshes' the build by deleting the necessary files
.PHONY: clean
clean:
	rm $(OBJS) memory ../memory












