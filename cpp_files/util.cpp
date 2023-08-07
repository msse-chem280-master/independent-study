#include <utility>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

#include "types.hpp"

// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

std::pair<Coordinates, double> read_xyz(std::string file_path)
{
    // Opens up a file stream forinput
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if (!infile.is_open()) {
        throw std::runtime_error("File path in read_xyz does not exist!");
    }

    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;

    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;

    // now the number of atoms
    infile >> num_atoms;

    // Uncomment to help troubleshoot
    // std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;

    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;

    for (int i = 0; i < num_atoms; i++) {
        AtomCoord coord;

        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];

        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}

double random_double(double lower_bound, double upper_bound)
{
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(re);
}

double random_integer(int lower_bound, int upper_bound)
{
    // dist will return [a,b] but we want [a,b)
    std::uniform_int_distribution<int> dist(lower_bound, upper_bound - 1);
    return dist(re);
}

void initialize_rng(void)
{
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
}