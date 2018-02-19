all: shared_lib


FC = gfortran
FFLAGS = -O3 -fPIC

CC = gcc
CFLAGS = -O3 -fPIC

#LIBS =  -llapack -lblas
#LIBS =  /usr/lib/liblapack.so.3 /usr/lib/libblas.so.3
#LIBDIR = /usr/lib/x86_64-linux-gnu/openblas/
LIBDIR = /usr/lib
LIBS =  $(LIBDIR)/liblapack.so.3 $(LIBDIR)/libblas.so.3

#OBJS      = fortran_stop_wrapper.o clock.o expokit.o mataid.o
OBJS      = fortran_stop_wrapper.o expokit.o 

.f.o:;  $(FC) $(FFLAGS) -c $<
.c.o:;  $(CC) $(CFLAGS) -c $<


# Figure out OS and architecture
OS = $(shell uname)
ifeq ($(OS), MINGW32_NT-6.1)
  OS=WINNT
endif

# file extensions
ifeq ($(OS), WINNT)
  SHLIB_EXT = dll 
else ifeq ($(OS), Darwin)
  SHLIB_EXT = dylib
else
  SHLIB_EXT = so
endif


shared_lib: libexpokit.$(SHLIB_EXT)

libexpokit.$(SHLIB_EXT): $(OBJS)
	$(FC) $(OBJS) -rdynamic -shared -o $@  $(LIBS)

clean:
	rm -f *.o *.$(SHLIB_EXT)

