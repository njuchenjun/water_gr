

ana.x:	test.o \
        shift_surf.o gwrite.o gwrite_poscar.o normal.o \
        pes_2023_c32.o pes_2023_int.o  \
        pes_h2o.o pes_eann.o energy1_au.o
	ifort -O2 $+ -o $@ -mkl

%.o:	%.f90
	ifort -O2 -c $< -o $@
%.o:    %.F90
	ifort -O2 -c $< -o $@
%.o:	%.f
	ifort -O2 -c $< -o $@


