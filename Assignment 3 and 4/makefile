# This is my make file for my ants program.
# I will call it 

# First, I will define the compiler to use.

CPP = g++ # We will use g++ to compile our files in the make file. 

# Compiler and linker flags.
CPPFLAGS = -pg -gstabs -O3 # Defines our compiler settings. We included optimization and the standard compiler to be used.
LDLIBS = -lboost_unit_test_framework # Defines our library flags. Needed for the make test/BOOST.

all: ants ants_original
# This line strings all of our libraries together (including the driver code) to form an executable, ants.
ants: ants.o single_time_step.o initial_speed.o initial_amount.o print.o ticktock.o
	${CPP} ${CPPFLAGS} -o ants ants.o single_time_step.o initial_speed.o ticktock.o \
	initial_amount.o print.o 

#The following codes create our library files from our source code.
#Note that the headers were not included - doing so will cause the make file to not compile properly!
ants.o: ants.cc
	${CPP} ${CPPFLAGS} -c -o ants.o ants.cc

initial_speed.o: initial_speed.cc
	${CPP} ${CPPFLAGS} -c -o initial_speed.o initial_speed.cc
	
initial_amount.o: initial_amount.cc
	${CPP} ${CPPFLAGS} -c -o initial_amount.o initial_amount.cc

simulation.o: simulation.cc
	${CPP} ${CPPFLAGS} -c -o single_time_step.o single_time_step.cc

print.o: print.cc
	${CPP} ${CPPFLAGS} -c -o print.o print.cc

ticktock.o: initial_speed.cc
	${CPP} ${CPPFLAGS} -c -o ticktock.o ticktock.cc

# Here is the part of the make file to compile the original ants file.
ants_original: ants_original.o ticktock.o
	${CPP} ${CPPFLAGS} -o ants_original ants_original.o ticktock.o 

ants_original.o: ants_original.cc
	${CPP} ${CPPFLAGS} -c -o ants_original.o ants_original.cc

# Lastly, we will include a clean command so that we can get rid of the files made in a previous run of the make file. This is executed by typing 'make clean' in the terminal.

clean: 
	rm -f ants.o ants_original.o initial_speed.o initial_amount.o single_time_step.o print.o ants ants ants_original ticktock.o gmon.out ants_test

# We have now included a test part of our make file. Running make test compiles the boost test code.

test:
	${CPP} ${CPPFLAGS} -o ants_test single_time_step.cc initial_speed.cc initial_amount.cc ants_test.cc ${LDLIBS}


