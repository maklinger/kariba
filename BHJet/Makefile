# Set the `GSL` variable to the base of your GSL installation
# If this happens to be already on a default path, no harm is done
# This option is mainly for semi-external installations, such as
# through Homebrew on macOS
GSL = /opt/homebrew/Cellar/gsl/2.8

# Path to Kariba root
KARIBA = ../Kariba

CXX = g++-15
CXXFLAGS = -std=c++17 -Wall -Wextra -fopenmp -fPIC -Ofast
CPPFLAGS = -I${GSL}/include -I${KARIBA}
LD = $(CXX)
LDFLAGS = -fPIC -L${GSL}/lib -fopenmp
LDLIBS = -lgsl -lm
RM = rm -f


SOURCES = bhjet.cpp EBL.cpp jetpars.cpp utils.cpp
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
