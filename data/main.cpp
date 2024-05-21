// main.cpp
#include "read_data.h"
#include "linalg.h"
#include "Markowitz.h"
#include <iostream>
using namespace std;
int  main (int  argc, char  *argv[])
{

    int numberAssets=83;
    int numberReturns=700;
    Lattice returnMatrix(numberAssets, Vector(numberReturns)); // a matrix to store the return data

    //read the data from the file and store it into the return matrix
    string fileName="asset_returns.csv";
    readData(returnMatrix,fileName);          // returnMatrix[i][j] stores the asset i, return j value

    Matrix daily_returns(returnMatrix);

    Lattice tr =  {
            {0.   , 0.005, 0.01 , 0.015, 0.02 , 0.025, 0.03 , 0.035, 0.04 ,
            0.045, 0.05 , 0.055, 0.06 , 0.065, 0.07 , 0.075, 0.08 , 0.085,
            0.09 , 0.095, 0.1}
    };

    Matrix target_returns(tr);

    Markowitz seven00(daily_returns, target_returns);
    Matrix results = seven00.results();
    results.shape();

    return 0;
}


/*// Example usage
int main() {
//
//    Matrix rowVector(Lattice(1, Vector(5, 0.0)));
//    Matrix columnVector(Lattice(5, Vector(1, 0.0)));
//    rowVector.prn();
//    rowVector.transpose().prn();

    // Define a simple symmetric positive definite matrix A
    Lattice A_data = {
            {2, 1},
            {1, 2}
    };
    Matrix A(A_data);

    // Define the vector b
    Lattice b_data = {
            {5},
            {6}
    };
    Matrix b(b_data);
    Lattice tmp =
    {
        {1.0, 0.0},
        {0.0, 1.0}
    };
    Matrix I(tmp);

    A.prn();
    b.prn();
    Matrix x  =  A.solver(b);
    x.prn();

    return 0;
}*/
