CONFIG = ../make.config
include $(CONFIG)

# Path to Kariba root
KARIBA = ../Kariba

ifndef GSL
	GSL_CPPFLAGS=
	GSL_LDFLAGS=
else
	GSL_CPPFLAGS=-I$(GSL)/include
	GSL_LDFLAGS=-L$(GSL)/lib
endif


CXXFLAGS = -std=c++17 -Wall -Wextra -fopenmp -fPIC -Ofast
CPPFLAGS = $(GSL_CPPFLAGS) -I${KARIBA}
LD = $(CXX)
LDFLAGS = -fPIC -fopenmp $(GSL_LDFLAGS)
LDLIBS = -lgsl -lm
RM = rm -f


SOURCES = bhjet.cpp jetpars.cpp utils.cpp
OBJECTS = $(subst .cpp,.o,$(SOURCES))
TARGET_OBJ = bhwrap.o
TARGET = bhwrap
PYOBJ = pyjetmain.o
PYLIB = pyjetmain.so
LIBKARIBA = ${KARIBA}/libkariba.a


.PHONY: all clean distclean


all: $(TARGET) $(PYLIB)

$(TARGET): $(OBJECTS) $(LIBKARIBA) $(TARGET_OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJECTS) $(TARGET_OBJ) $(LIBKARIBA) $(LDLIBS)

$(PYLIB): $(OBJECTS) $(LIKARIBA) $(PYOBJ)
	$(LD) $(LDFLAGS) -shared -o $(PYLIB) $(OBJECTS) $(PYOBJ) $(LIBKARIBA) $(LDLIBS)

clean:
	$(RM) $(OBJECTS) $(KARIBA_OBJECTS) $(TARGET_OBJ) $(PYOBJ)

distclean: clean
	$(RM) $(TARGET) $(PYLIB)
