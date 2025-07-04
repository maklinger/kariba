CONFIG = ../make.config
include $(CONFIG)

KARIBA = ../Kariba


ifndef GSL
	GSL_CPPFLAGS=
	GSL_LDFLAGS=
else
	GSL_CPPFLAGS=-I$(GSL)/include
	GSL_LDFLAGS=-L$(GSL)/lib
endif


# Path to Kariba root
KARIBA = ../Kariba

CXXFLAGS = -std=c++17 -Wall -Wextra -fopenmp -fPIC -Ofast
CPPFLAGS = $(GSL_CPPFLAGS) -I${KARIBA}
LD = $(CXX)
LDFLAGS = -fPIC -fopenmp $(GSL_LDFLAGS)
LDLIBS = -lgsl -lm
RM = rm -f


SOURCES = utils.cpp
OBJECTS = $(subst .cpp,.o,$(SOURCES))

TARGET_CORONA = kariba_corona_examples
TARGET_PARTICLES = kariba_particles_examples
TARGET_SINGLEZONE = kariba_singlezone_examples
OBJECT_CORONA = kariba_corona_examples.o
OBJECT_PARTICLES = kariba_particles_examples.o
OBJECT_SINGLEZONE = kariba_singlezone_examples.o

LIBKARIBA = ${KARIBA}/libkariba.a


.PHONY: all clean distclean


all: $(TARGET_CORONA) $(TARGET_PARTICLES) $(TARGET_SINGLEZONE)


$(TARGET_CORONA): $(OBJECTS) $(LIBKARIBA) $(OBJECT_CORONA)
	$(LD) $(LDFLAGS) -o $(TARGET_CORONA) $(OBJECTS) $(OBJECT_CORONA) $(LIBKARIBA) $(LDLIBS)

$(TARGET_PARTICLES): $(OBJECTS) $(LIBKARIBA) $(OBJECT_PARTICLES)
	$(LD) $(LDFLAGS) -o $(TARGET_PARTICLES) $(OBJECTS) $(OBJECT_PARTICLES) $(LIBKARIBA) $(LDLIBS)

$(TARGET_SINGLEZONE): $(OBJECTS) $(LIBKARIBA) $(OBJECT_SINGLEZONE)
	$(LD) $(LDFLAGS) -o $(TARGET_SINGLEZONE) $(OBJECTS) $(OBJECT_SINGLEZONE) $(LIBKARIBA) $(LDLIBS)


clean:
	$(RM) $(OBJECTS) $(OBJECT_CORONA) $(OBJECT_PARTICLES) $(OBJECT_SINGLEZONE)

distclean: clean
	$(RM) $(TARGET_CORONA) $(TARGET_PARTICLES) $(TARGET_SINGLEZONE)
