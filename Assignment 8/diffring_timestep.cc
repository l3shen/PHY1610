// Time-step module for diffusion model.

#include "diffring_timestep.h"
#include <cblas.h>

using namespace std;

// The following two functions were blank, and were needed to be filled in. 

// This function is used to perform a single time step in the density field.
// The F matrix is the derivative matrix that describes time evolution, and the P matrix is the density matrix itself.
void perform_time_step(const rarray<double,2>& F, rarray<double,1>& P) {

	// Initialize an array to store the calculation during CBLAS.
	rarray<double,1> P_new(P.extent(0)); 

	// rarray method to find the number of elements in in each dimension. Needed for CBLAS; gives order.
	int n = F.extent(0); // Could not use 'size.()' as it lead to a segmentation fault.

	// Initialize scalar parameters for our cblas call.
  	float alpha = 1.0; 
  	float beta = 0.0; 

	// Use cblas_dsymv to compute a matrix-vector product for a symmetric matrix (which ours conveniently is).
  	cblas_dsymv(CblasRowMajor, CblasUpper, n, alpha, F.data(), n, P.data(), 1, beta, P_new.data(), 1);

	// Place the results of our calculation in to our probability matrix, passed by reference.
	P = P_new.copy();

}

// This fills the F matrix with its appropriate items.
void fill_time_step_matrix(rarray<double,2>& F, double D, double dt, double dx) {

	// Similarly to above, needed for CBLAS, determines order.
  	int o = F.extent(0);
  
	// Fill our matrix in order to reserve memory (?) - gets Segmentation Fault without.
  	F.fill(0.0);

	// Calculating the elements of our matrix to then fill.
  	double offdiagonal = (D*dt)/(pow(dx,2)); 
  	double diagonal = 1-(2*offdiagonal);     

	// Initiating a loop to fill our rarray.
  	for (int i = 1; i <= o - 2; i++) {

		// Fill our diagnonal elements.
		F[i][i] = diagonal; 

		// Fill our offdiagonal elements.
    		F[i][i - 1] = offdiagonal;  
    		F[i][i + 1] = offdiagonal;

  	}
	
	// Filling our matrix.
	// Offdiagonal elements.
  	F[0][o - 1] = offdiagonal;
  	F[o - 1][0] = offdiagonal;
  	F[o - 1][o - 2] = offdiagonal;
  	F[o - 2][o - 1] = offdiagonal;
  	F[0][1] = offdiagonal;

	// Diagonal elements.
  	F[0][0] = diagonal;
  	F[o - 1][o - 1] = diagonal;  	
}
