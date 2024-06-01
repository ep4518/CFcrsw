//write_data.h
#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#include <linalg.h>
#include <fstream>
#include <sstream>
#include <unistd.h>

struct Result {
    int id;
    Matrix back_test;
    Matrix weights;
};

// Write results array to a csv format
void write_data(const Result *results, size_t size);

// converts a Lattice to a string for write_data usage
std::string vectorToString(const Lattice& M);

#endif