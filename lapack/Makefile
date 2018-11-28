FC = gfortran

FFLAGS = -fdefault-integer-8 -I${F95ROOT}/include/intel64/ilp64 -m64 -I${MKLROOT}/include

LDFLAGS = ${F95ROOT}/lib/intel64/libmkl_blas95_ilp64.a ${F95ROOT}/lib/intel64/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

DEPS =

OBJ = lapack_matrix_inversion.o

%.o: %.f08 $(DEPS)
	$(FC) -c -o $@ $< $(FFLAGS) 

a.out: $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)

.PHONY: clean, veryclean

clean:
	rm -f *.o

veryclean:
	rm -f *.o a.out

# $@ means left of :
# $^ means right of :
# $< means the firts item on the right of :
# do not look for a filenames on the right of .PHONY:
