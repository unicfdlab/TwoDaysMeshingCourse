LAMATHDIR = $(LaMathHome)
CFLAGS = -Wall -O3  -I ../../LaMath/SPOOLES -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE
FFLAGS = -Wall -O3

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.6.1.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

SPOOLES_DIR=$(LAMATHDIR)/SPOOLES

LIBS = \
       $(SPOOLES_DIR)/spooles.a \
	$(LAMATHDIR)/ARPACK/libarpack_linux64.a \
       -lpthread -lm -lc

ccx_2.6.1: $(OCCXMAIN) ccx_2.6.1.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.6.1.c; $(FC) -Wall -O3 -o $@ $(OCCXMAIN) ccx_2.6.1.a $(LIBS)

ccx_2.6.1.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

