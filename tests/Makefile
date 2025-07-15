CONFIG = ../make.config
include $(CONFIG)

# Path to Kariba root
KARIBA = ../src


ifndef GSL
	GSL_CPPFLAGS=
	GSL_LDFLAGS=
else
	GSL_CPPFLAGS=-I$(GSL)/include
	GSL_LDFLAGS=-L$(GSL)/lib
endif


CXXFLAGS += $(OPENMP)
CPPFLAGS = -I${GSL}/include -I${KARIBA} -I.
LD = $(CXX)
LDFLAGS = -L${GSL}/lib $(OPENMP) $(EXTRA_LDFLAGS)
LDLIBS = -lgsl -lm
RM = rm -f

# Kariba library
LIBKARIBA = ${KARIBA}/libkariba.a


SOURCES = test_bknpower.cpp test_compton.cpp test_cyclosyn.cpp test_distributions.cpp test_ebl.cpp test_particles.cpp test_powerlaw.cpp test_radiation.cpp
OBJECTS = $(subst .cpp,.o,$(SOURCES))
MAIN_OBJ = test_main.cpp
TEST_MAIN = test_main


.PHONY: all tests main clean distclean

all: main

main: $(OBJECTS) $(LIBKARIBA) $(MAIN_OBJ)
	$(LD) $(LDFLAGS) -o $(TEST_MAIN) $(OBJECTS) $(MAIN_OBJ) $(LIBKARIBA) $(LDLIBS)


clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(ALL_TESTS)
