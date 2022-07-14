CXX = mpicxx -std=c++11

# production run
CXXFLAGS = -Wall -O3 -ffast-math -ffinite-math-only -fno-trapping-math -fpermissive -march=corei7 -DVIS_PDF
#CXXFLAGS = -cxx=icpc -Wall -ipo -O3 -no-prec-div -fp-model fast=2 -xSSE4.2



# for debugging
#CXXFLAGS = -Wall -g -O0 -fpermissive

PROG0 = hybrid_sim

INCLUDES = -I.

LIBS = -L.

OBJECTS0 = geo_coord.o \
           hybrid_particles.o \
           hybrid_distr.o \
           block_thomas.o \
           hybrid_fields.o \
           hybrid_moments.o \
           $(PROG0).o

all: $(PROG0)

$(PROG0): $(OBJECTS0)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(PROG0) $(OBJECTS0) $(LIBS)

clean:
	rm -rf *.o $(PROG0) core

%.o: %.cc *.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -c $<

