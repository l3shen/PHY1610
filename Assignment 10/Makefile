#
# Makefile for 1dwave code for SciNet's High Performance Scientific Computing
#  
# Modules required
#
#module load gcc/4.8.1
#module load pgplot/5.2.2-gcc


CXX=mpic++
CPPFLAGS=-I${SCINET_PGPLOT_INC}
CXXFLAGS=-g -O0 -std=c++11 -fopenmp
LDLIBS=-lcpgplot -lpgplot -lX11 -lgfortran -lpng -fopenmp 
LDFLAGS=-L${SCINET_PGPLOT_LIB}

all: wave1d

# Makefile for wave equation

wave1d: wave1d.o ticktock.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

wave1d_pglot.o: wave1d.cc ticktock.h
	$(CXX) -c -fopenmp $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

ticktock.o: ticktock.cc ticktock.h
	$(CXX) -c -fopenmp $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

wave1d_mpi: wave1d_mpi_final.o ticktock.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

wave1d_mpi_no_out: wave1d_mpi_final_nooutput.o ticktock.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.o
