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

void write_data(const Result *results, size_t size);

std::string vectorToString(const Lattice& M);

#endif