# Makefile for the signal analysis program.
# Run make && make install, followed by ./signal_analysis. 

CPP = g++ # We will use g++ to compile our files in the make file. 

# Compiler and linker flags.
CPPFLAGS = -I${BLAS_INC} -O2 # Using O2, albeit O3 would probably be alright.
LDLIBS = -lfftw3 -L${BLAS_LIB} -lopenblas # Defines our library flags. Needed for FFTW3 and OpenBLAS.

all: signal_analysis

signal_analysis: driver.o fftw.o filename.o file_reader.o power_spectrum_cc.o highest_value.o
	${CPP} ${CPPFLAGS} -o signal_analysis driver.o fftw.o filename.o file_reader.o \
	power_spectrum_cc.o highest_value.o ${LDLIBS}

driver.o: driver.cc
	${CPP} ${CPPFLAGS} -c -o driver.o driver.cc ${LDLIBS}

fftw.o: fftw.cc
	${CPP} ${CPPFLAGS} -c -o fftw.o fftw.cc ${LDLIBS}
	
filename.o: filename.cc
	${CPP} ${CPPFLAGS} -c -o filename.o filename.cc ${LDLIBS}

file_reader.o: file_reader.cc
	${CPP} ${CPPFLAGS} -c -o file_reader.o file_reader.cc ${LDLIBS}

power_spectrum_cc.o: power_spectrum_cc.cc
	${CPP} ${CPPFLAGS} -c -o power_spectrum_cc.o power_spectrum_cc.cc ${LDLIBS}

highest_value.o: highest_value.cc
	${CPP} ${CPPFLAGS} -c -o highest_value.o highest_value.cc ${LDLIBS}

clean: 
	rm -f driver.o fftw.o filename.o file_reader.o power_spectrum_cc.o highest_value.o signal_analysis


