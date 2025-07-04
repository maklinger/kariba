CONFIG = ../make.config
include $(CONFIG)


ifndef GSL
	GSL_CPPFLAGS=
	GSL_LDFLAGS=
else
	GSL_CPPFLAGS=-I$(GSL)/include
	GSL_LDFLAGS=-L$(GSL)/lib
endif


CXXFLAGS = -std=c++17 -Wall -Wextra -fPIC -O3 -fopenmp
CPPFLAGS = $(GSL_CPPFLAGS)
LD = $(CXX)
LDFLAGS = -fPIC -fopenmp $(GSL_LDFLAGS)
LDLIBS = -lgsl
AR = ar
AR_OPTS = rcs
RM = rm -f


SOURCES = BBody.cpp Bknpower.cpp Compton.cpp Cyclosyn.cpp EBL.cpp Electrons.cpp GammaRays.cpp Kappa.cpp Mixed.cpp Neutrinos_pg.cpp Neutrinos_pp.cpp Particles.cpp Powerlaw.cpp Radiation.cpp ShSDisk.cpp Thermal.cpp
OBJECTS = $(subst .cpp,.o,$(SOURCES))
LIBSHARED = libkariba.so
LIBSTATIC = libkariba.a


.PHONY: all clean distclean


all: shared static

# Use .so on macOS as well; the difference between .so and .dylib is minimal
# This allows the use of the same extension and same -shared option
shared: $(OBJECTS)
	$(LD) $(LDFLAGS) -shared -o $(LIBSHARED) $(OBJECTS) $(LDLIBS)

# macOS compatible way to create a static library
static: $(OBJECTS)
	$(AR) $(AR_OPTS) $(LIBSTATIC) $(OBJECTS)

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(LIBSHARED) $(LIBSTATIC)
