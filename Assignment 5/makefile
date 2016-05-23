# This is my make file for my minima finding program.

# First, I will define the compiler to use.

CPP = g++ # We will use g++ to compile our files in the make file. 

# Compiler and linker flags. For both local and SciNet compilation.
CPPFLAGS = -I${GSLINC} -O2
LDLIBS = -lgsl -lgslcblas
SCINET_CPPFLAGS = -I${SCINET_GSL_INC} -O2
SCINET_LDFLAGS = -L${SCINET_GSL_LIB}

# For local machine; just run make.
all: minima_and_mass_ext

minima_and_mass_ext: main_driver.o pot_energy.o f_all_min.o 
	${CPP} -o minima_and_mass_ext main_driver.o pot_energy.o f_all_min.o ${LDLIBS}

main_driver.o: main_driver.cc
	${CPP} ${CPPFLAGS} -c -o main_driver.o main_driver.cc

pot_energy.o: pot_energy.cc
	${CPP} ${CPPFLAGS} -c -o pot_energy.o pot_energy.cc
	
f_all_min.o: f_all_min.cc
	${CPP} ${CPPFLAGS} -c -o f_all_min.o f_all_min.cc


# For SciNet, run make scinet.
scinet: minima_and_mass_ext_scinet

minima_and_mass_ext_scinet: main_driver.o_scinet pot_energy.o_scinet f_all_min.o_scinet 
	${CPP} ${SCINET_CPPFLAGS} -o minima_and_mass_ext main_driver.o pot_energy.o f_all_min.o ${SCINET_LDFLAGS} ${LDLIBS}

main_driver.o_scinet: main_driver.cc
	${CPP} ${SCINET_CPPFLAGS} -c -o main_driver.o main_driver.cc

pot_energy.o_scinet: pot_energy.cc
	${CPP} ${SCINET_CPPFLAGS} -c -o pot_energy.o pot_energy.cc
	
f_all_min.o_scinet: f_all_min.cc
	${CPP} ${SCINET_CPPFLAGS} -c -o f_all_min.o f_all_min.cc

# To clean files.

clean: 
	rm -f minima_and_mass_ext main_driver.o pot_energy.o f_all_min.o global_minima.txt


