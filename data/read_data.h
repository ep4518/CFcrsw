//read_data.h
#ifndef _READ_DATA_H
#define _READ_DATA_H
#include <sstream>
#include <vector>

double string_to_double( const std::string& s );
void readData(std::vector<std::vector<double> >& returnMatrix,std::string fileName);

#endif