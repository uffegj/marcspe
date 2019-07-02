#F90          = /software/fortran/intel/fce/9.0/bin/ifortbin
#F90LINKER    = /software/fortran/intel/fce/9.0/bin/ifortbin
F90          = ifort
F90LINKER    = ifort
#F90         = f95
#F90LINKER   = f95
#F90         = gfortran
#F90LINKER   = gfortran

### options ###

# for IFORT-Compiler
#===================
DEFS      = -D IFORT
#FFLAGS    = -r8 -C -fpe0 -g -traceback
FFLAGS    = -r8 -O3 -g -traceback
INCLUDES  = 
LFLAGS    = $(FFLAGS)
# -L/home/woitke/data/LAPACK/mkl_9.0_beta/lib/32 -lm -lmkl -lmkl_lapack 
# -L /disks/theory1/woitke/LAPACK3.0 -llapack 

# for NAG-Compiler (geht nicht auf delphinus)
#============================================
#DEFS      = 
#FFLAGS    = -fixed -r8 -g -w=all
#INCLUDES  = 
#LFLAGS    = -fixed -r8 -g

# for gfortran-Compiler
#======================
#DEFS      = 
#FFLAGS    = -fdefault-real-8 -g -O5
#INCLUDES  = 
#LFLAGS    = -fdefault-real-8 -g -O5

EXECS     = ../static_weather8

OBJECTS  = drift_data.o \
           INIT_AMES.o \
           INIT_DUSTCHE.o \
           INIT_PHOENIX.o \
           INIT_TestCases.o \
           INIT_MARCS.o \
           INIT.o \
           THERMO.o \
           SUPERSAT.o \
           EFF_SUPERSAT.o \
           NUCLEATION.o \
           KLASS_NUK_TI.o \
           KLASS_NUK_C.o \
           JEVAP.o \
           CHI_NET.o \
           ABSCHLUSS_BE.o \
           dist_from_moments.o \
           static_weather8.o \
           FF.o \
           STATIC_EPS.o \
           ELFUNK.o \
           limex.o \
           limex_dense.o \
           limex_aux.o \
           FCN.o \
           FDJAC.o \
           FUNCV.o \
           JACOBI.o \
           GAUSS.o \
           GGCHEM.o \
           smchem.o \
           DISPOL.o \
           stindex.o \
           KOPF.o \
           OUTPUT.o \
           is_nan.o \
           slatec_routines.o \
	   DESTIFF.o

default: sw

#
sw:      $(OBJECTS)
	 $(F90LINKER) -o $(EXECS) $(OBJECTS) $(LFLAGS) 
#
clean:
	-/bin/rm -f *.o 

#
realclean:
	-/bin/rm -f *.o *~ *.mod $(EXECS)

.f.o:
	$(F90) $(FFLAGS) -c $<
.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

