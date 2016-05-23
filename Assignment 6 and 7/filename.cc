#include "filename.h"

// This function generates a string to use as a filename when reading the data.

const std::string Filename(int a) {
		std::string filename; // Our final returned string.
		const std::string filename_1 = "detection"; // The common first part.
		const std::string filename_3 = ".rat"; // The common last part.

		// Generate a string from an integer.
		int value = a;
		std::string filename_2;
		std::ostringstream convert;
		convert << value;
		filename_2 = convert.str();

		// Since the naming convention for files that are <10 is 0x, we will add a zero to those cases to match up with the file names.
		if (a < 10) {
			std::string zero = "0";
			filename_2 = zero + filename_2;
		}

		// Return our final file name.
		return filename = filename_1 + filename_2 + filename_3;
}
