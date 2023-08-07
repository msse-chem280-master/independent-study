#pragma once

#include <string>

#include "types.hpp"

/*! Reads in an XYZ file from the given path */
std::pair<Coordinates, double> read_xyz(std::string file_path);


/*! Generate a random double within a given range */
double random_double(double lower_bound, double upper_bound);


/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double random_integer(int lower_bound, int upper_bound);


/*! Initialize the random number generator */
void initialize_rng(void);