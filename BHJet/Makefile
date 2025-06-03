# Set the `GSL` variable to the base of your GSL installation
# If this happens to be already on a default path, no harm is done
# This option is mainly for semi-external installations, such as
# through Homebrew on macOS
GSL = /opt/homebrew/Cellar/gsl/2.8

CXX = g++-15
CXXFLAGS = -std=c++17 -fopenmp -fPIC -Ofast
CPPFLAGS = -I${GSL}/include
LD = $(CXX)
LDFLAGS = -fPIC -L${GSL}/lib -fopenmp
LDLIBS = -lgsl -lm
RM = rm -f


SOURCES = bhjet.cpp bhwrap.cpp EBL.cpp jetpars.cpp utils.cpp
OBJECTS = $(subst .cpp,.o,$(SOURCES))
TARGET = bhjet
LIBKARIBA = ../Kariba/libkariba.a


.PHONY: all clean distclean


all: $(TARGET)

$(TARGET): $(OBJECTS) $(LIBKARIBA)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJECTS) $(LIBKARIBA) $(LDLIBS)

clean:
	$(RM) $(OBJECTS) $(KARIBA_OBJECTS)

distclean: clean
	$(RM) $(TARGET)
