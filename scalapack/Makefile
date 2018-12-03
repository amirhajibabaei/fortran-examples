FC = mpifort

FFLAGS = -fdefault-integer-8 -I${F95ROOT}/include/intel64/ilp64 -m64 -I${MKLROOT}/include -Wall -Wextra -ffree-line-length-200 

LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_ilp64 -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_ilp64 -lpthread -lm -ldl

DEPS =

OBJ = blacs_scalapack_basics.o 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  rules
a.out: $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)


%.o: %.f08 $(DEPS)
	$(FC) -c -o $@ $< $(FFLAGS) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.PHONY: clean, veryclean

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod a.out


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ notes
# $@ means left of :
# $^ means right of :
# $< means the firts item on the right of :
# do not look for a filenames on the right of .PHONY:
