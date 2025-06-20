# Set the `GSL` variable to the base of your GSL installation
# If this happens to be already on a default path, no harm is done
# This option is mainly for semi-external installations, such as
# through Homebrew on macOS
GSL = /opt/homebrew/Cellar/gsl/2.8

CXX = g++-15
CXXFLAGS = -std=c++17 -Wall -Wextra -fPIC -Ofast
CPPFLAGS = -I${GSL}/include -I.
LD = $(CXX)
LDFLAGS = -fPIC -L${GSL}/lib
LDLIBS = -lgsl
AR = ar
AR_OPTS = rcs
RM = rm -f


SOURCES = BBody.cpp Bknpower.cpp Compton.cpp Cyclosyn.cpp Electrons.cpp GammaRays.cpp Kappa.cpp Mixed.cpp Particles.cpp Powerlaw.cpp Radiation.cpp ShSDisk.cpp Thermal.cpp
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
