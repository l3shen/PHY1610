//
// wave1d.cc - Simulates a one-dimensional damped wave equation modified
// to use pgplot for runtime display of solution with 1 second between frames.  
// 
// SciNet - March 2015

// Modified for Assignment 10.

#include <iostream>
#include <ostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <rarray>
#include <cpgplot.h>
#include "ticktock.h"
#include "inifile.h"
#include <mpi.h>
#include <string>

int main(int argc, char* argv[])
{

    // Initializing MPI variables.
    int ierr;
    int rank, size;

    // Open inifile and parse (using Inifile class from inifile.h)
    Inifile parameter(argc==1?"default.txt":argv[1]);

    // Physical parameters
    double  c       = parameter.get<double>("c", 1.0);     // wave speed
    double  tau     = parameter.get<double>("tau", 20.0);  // damping time
    double  x1      = parameter.get<double>("x1", -26.0);  // left most x value
    double  x2      = parameter.get<double>("x2", +26.0);  // right most x value

    // Simulation parameters
    double  runtime = parameter.get<double>("runtime", 50.0);   // how long should the simulation try to compute?
    double  dx      = parameter.get<double>("dx", 0.01);        // spatial grid size  

    // Output parameters
    double  outtime =  parameter.get<double>("outtime", 1.0); // how often should a snapshot of the wave be written out? 

    bool    graphics = parameter.get<bool>("graphics", true);   // output to graphics (with 1 sec delay)  or to a file?

    // Derived parameters
    int     ngrid   = (x2-x1)/dx - 1;  // number of x points, modified due to potential error?
    // This was modified from the original!
    int     npnts_global   = ngrid + 2;   // number of x points including boundary points
    double  dt      = 0.5*dx/c;    // time step size
    int     nsteps  = runtime/dt;  // number of steps of that size to reach runtime
    int     nper    = outtime/dt;  // how many step s between snapshots

    // Report all the values.
    std::cout << "#c        " << c       << std::endl;
    std::cout << "#tau      " << tau     << std::endl;
    std::cout << "#x1       " << x1      << std::endl;
    std::cout << "#x2       " << x2      << std::endl;
    std::cout << "#runtime  " << runtime << std::endl;
    std::cout << "#dx       " << dx      << std::endl;
    std::cout << "#outtime  " << outtime << std::endl; 
    std::cout << "#ngrid    " << ngrid   << std::endl;
    std::cout << "#dt       " << dt      << std::endl;
    std::cout << "#nsteps   " << nsteps  << std::endl;    
    std::cout << "#nper     " << nper    << std::endl;
    std::cout << "#graphics " << int(graphics) << std::endl;

    // Initializing MPI.
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculating size of local domain, and number of points per local grid.
    int local_ngrid = 0;
    int modulus = ngrid % size;

	// Accounting for the modulus, such that certain processors carry more digits (i.e. not all local grids are same size).
    if (modulus == 0) {
	local_ngrid = (ngrid)/size;
    }
    // Accounting for local cells that will carry an extra grid point.
    else if (modulus !=0) {
        local_ngrid = (ngrid - modulus)/size;
	if ((rank >= 0) and (rank <= (modulus - 1)))
	    local_ngrid += 1;
    }

    // Determining the local domain by including the ghost cells (accounting for the boundaries).
    int local_domain = local_ngrid + 2;

    // Reassigning local domain to npnts to avoid having to rewrite everything downstream.
    int npnts = local_domain;

    // Define and allocate arrays.
    rarray<float,1> rho_prev(npnts); // time step t-1
    rarray<float,1> rho(npnts);      // time step t
    rarray<float,1> rho_next(npnts); // time step t+1
    rarray<float,1> rho_init(npnts); // initial values
    rarray<float,1> x(npnts);        // x values

    // Initialize (MPI).
	// Find local edges (not boundaries).
	int x1_val = rank * (local_ngrid)+1;
	int x2_val = x1_val + (local_ngrid-1);
	// Ensure that edge is not part of the larger ngrids brought on by a non-zero remainder.
	if (rank >= modulus) {
		x1_val += modulus;
		x2_val += modulus;
	}

	// Determine the local boundaries for the local grid.
	double local_x1 = x1 + (x1_val)*dx;
	double local_x2 = x1 + (x2_val)*dx;

	// Determine the x array (i.e. positions of the grid).
	for (int i = 0; i <= (npnts - 1); i++) {
	  x[i] = local_x1 + (i-1)*dx;
	}

	// Excite; generate zero rho grid first.
	for (int i = 0 ; i < npnts; i++) {
        rho[i] = 0.0;
        rho_prev[i] = 0.0;
        rho_next[i] = 0.0;
    }

	// Excite (MPI) - fill in the non-zero elements according to the local grid size.
    for (int i = 0; i < npnts ; i++){
		if(((i+x1_val-1) >= ((npnts_global/4)+1)) and ((i+x1_val-1) < (3*npnts_global/4))){
			rho[i] = 0.25 - fabs(float((i+x1_val-1)-npnts_global/2)/float(npnts_global));
        	rho_prev[i] = rho[i];
        	rho_init[i] = rho[i];
		}
    } 

    // Plot or Write out data.
    std::ofstream dataFile;
    int red, grey, white;

    if (graphics) {
       cpgbeg(0, "/xwindow", 1, 1);
       cpgask(0);
       red = 2; cpgscr(red,1.,0.,0.);
       grey = 3; cpgscr(grey,.2,.2,.2);
       white = 4; cpgscr(white,1.0,1.0,1.0);
       cpgsls(1); cpgslw(6); cpgsci(white);
       cpgslw(2);
       cpgenv(x1, x2, 0., 0.25, 0, 0);
       cpglab("x", "rho", "Wave Test");
       cpgsls(1); cpgslw(6); cpgsci(white);
       cpgline(npnts, x.data(), rho.data());
       cpgsls(2); cpgslw(12); cpgsci(red);
       cpgline(npnts, x.data(), &rho_init[0]);
    } else {
	   // Determine output file name (unique for each rank).
       std::string fn_one = "dataFilename";
	   std::string fn_three = ".out";
       std::ostringstream convert;
       convert << rank;
       std::string fn_two = convert.str();
       const std::string dataFilename = fn_one + fn_two + fn_three;
	   convert.clear();

       dataFile.open(dataFilename.c_str());
       dataFile << nper << ","   
                << npnts       << "\n";
       dataFile << time << "\n";

       // Take care of guard cells.
       int guard_1 = 1;
       int guard_2 = (npnts-1);

	   // Dispose of guard cells when necessary.
       if (rank == 0) {
       		guard_1 = 0;
       }
	   else if (rank == (size - 1)) {
            guard_2 = (npnts);
       }

       for (int i = guard_1; i < guard_2; i++ ) 
          dataFile << x[i] << " " << rho[i] << " \n";  
       dataFile << "\n";
    }

    // measure time
    TickTock tt;
    tt.tick();

    // Initialize communication.
	MPI_Request request[4];
	MPI_Status status[4];

	int left = rank - 1;
	int right = rank + 1;

	// Ensure boundary coniditions for communication.
    if (right >= size ) {
        right = MPI_PROC_NULL;
    }
	if (left < 0) {
		left = MPI_PROC_NULL;
    }


	// Take timesteps (MPI)
    for (int s = 0; s < nsteps; s++) {

		// Communication, a la lecture 2 of MPI.
		if(s!=0){
		// Send right, receive left.
		  ierr = MPI_Isend(&rho[local_ngrid], 1, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &request[0]);
          ierr = MPI_Irecv(&rho[0], 1, MPI_FLOAT, left, 1, MPI_COMM_WORLD, &request[1]);
        // Send left, receive right.
          ierr = MPI_Isend(&rho[1], 1, MPI_FLOAT, left, 2, MPI_COMM_WORLD, &request[2]);
          ierr = MPI_Irecv(&rho[local_ngrid + 1], 1, MPI_FLOAT, right, 2, MPI_COMM_WORLD, &request[3]);
		}

        // Set zero dirichlet boundary conditions, modified for MPI.
		if (rank == 0) {
           rho[0] = 0.0;
		}
        else if (rank == (size - 1)) {
           rho[npnts-1] = 0.0;
		}

		// Do not communicate while calculating.
		if (s!=0) {
			ierr = MPI_Waitall(4, request, status);
		}

        	// Evolve (MPI)
		// Calculate boundary values.
		// First boundary.
		int i=1;
        float laplacian = pow(c/dx,2)*(rho[i+1] + rho[i-1] - 2*rho[i]);
        float friction = (rho[i] - rho_prev[i])/tau;
        rho_next[i] = 2*rho[i] - rho_prev[i] + dt*(laplacian*dt-friction);

		// Last boundary.
        int z=local_ngrid;
        laplacian = pow(c/dx,2)*(rho[z+1] + rho[z-1] - 2*rho[z]);
        friction = (rho[z] - rho_prev[z])/tau;
        rho_next[z] = 2*rho[z] - rho_prev[z] + dt*(laplacian*dt-friction);

		// Calculate non-boundary values.
        for (int i = 2; i <= (local_ngrid - 1); i++) {
            float laplacian = pow(c/dx,2)*(rho[i+1] + rho[i-1] - 2*rho[i]);
            float friction = (rho[i] - rho_prev[i])/tau;
            rho_next[i] = 2*rho[i] - rho_prev[i] + dt*(laplacian*dt-friction);
        }

        // Rotate array pointers so t+1 becomes the new t etc.
        rarray<float,1> temp;
        temp     = rho_prev;
        rho_prev = rho;
        rho      = rho_next;
        rho_next = temp;


        //Output every nper
        if ((s+1)%nper == 0) {
           if (graphics) {
              cpgbbuf();
              cpgeras();
              cpgsls(1); cpgslw(6); cpgsci(white);
              cpgslw(2);
              cpgenv(x1, x2, 0., 0.25, 0, 0);
              cpglab("x", "rho", "Wave test");  //t=s*dt
              cpgsls(2); cpgslw(12); cpgsci(red);
              cpgline(npnts, x.data(), rho.data());
              cpgsls(1); cpgslw(6); cpgsci(white);
              cpgline(npnts, x.data(), rho_init.data());
              cpgebuf();
              sleep(1); // artificial delay!
           } else {
              	dataFile << time << "\n";
				// Avoid outputting guard cells! Similar to first instance of file sending.
				int guard_1 = 1;
				int guard_2 = (npnts-1);

				if( rank == 0){
					guard_1 = 0;
				} 
				else if (rank == (size - 1)) {
					guard_2 = npnts;
				}
				// Same setup as before.
              	for (int i = guard_1; i < guard_2; i++)
                 	dataFile<< x[i] << " " << rho[i] << "\n";
              	dataFile << "\n";
           }
        }

    }

    // Output measured runtime.
    std::cout << "Walltime = " << tt.silent_tock() << " sec."  << std::endl;

    // Close file.
    if (not graphics)
       dataFile.close();

	// Close MPI.
   	ierr = MPI_Finalize();

   	return 0;
}
