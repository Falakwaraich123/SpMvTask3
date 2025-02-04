TARGETS = gsl mkl
CC = gcc
ICC = icx
OFLAGS = -g3 -O3 -march=native -mtune=icelake-server
WARNFLAGS = -Wall -Wextra
IVEC = -vec -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize 
CVEC = -fopt-info-vec-optimized -ffast-math 
CFLAGS =  $(OFLAGS) $(WARNFLAGS) $(CVEC) -D_GSL_
IFLAGS =  $(OFLAGS) $(WARNFLAGS) $(IVEC) -D_MKL_ 
LDLIBS = -lm -lgsl -lgslcblas
ILDLIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl

SRCS = timer.c my_dense.c my_csr.c my_coo.c my_csc.c spmv.c
IOBJS = $(SRCS:.c=.o)
OBJS = $(SRCS:.c=_g.o)

all: $(TARGETS)

gsl: $(OBJS) 	
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@
mkl: $(IOBJS)
	$(ICC) $(IFLAGS) $^ $(ILDLIBS) -o $@

$(IOBJS): $(SRCS)
	$(ICC) $(IFLAGS) $^ -c

%_g.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

cleanall:
	$(RM) $(OBJS) $(IOBJS) $(TARGETS) *~

