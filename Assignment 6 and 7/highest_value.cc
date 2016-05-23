#include "highest_value.h"

// This function returns the N_Highest number of highest terms in an array.
rarray<double,1> Highest_Value(rarray<double,1> Coefficients, rarray<int,1> &File_Index, int N_Highest) {
	
	rarray<double,1> CurrentHighVal(N_Highest); // Generate an array to fill.
	CurrentHighVal.fill(0.0); // Begin by filling it with zero; as we find greater values, we will change this.
	double Limit = 1.1; // This is our limit for searching for a maximum; 1 is the highest maximum, so we assure the first value cannot be greater than 1.

	// Define a for loop to find the 5 highest values, but can be changed to x highest values instead.
	int start = 0;
	int end = (Coefficients.size() - 1);
	for (int j = start; j <= (N_Highest - 1); j++) {
		// This for loop will find the highest value in the array and append that value to the i-th index of the array.
		for (int i = start; i <= end; i++) {
			if ((Coefficients[i] > CurrentHighVal[j]) and (Coefficients[i] < Limit)) {
				CurrentHighVal[j] = Coefficients[i];
				File_Index[j] = i;
			}
		}

		// Change the limit such that the next search finds the next highest value.
		Limit = CurrentHighVal[j];
	}

	// Return something we can work with.
	return CurrentHighVal;

}

// This is our print function. Nothing too complex here, so I won't really explain it.
void Print(rarray<int,1> File_Index, rarray<double,1> HighestValues, int N_Max) {
	for (int a = 0; a <= (N_Max - 1); a++) {
		std::cout << "#" << a+1 << "; The file number is: " << File_Index[a] << " The correlation coefficient is: " << HighestValues[a] << std::endl;
	}
}
