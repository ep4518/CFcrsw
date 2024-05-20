#include "read_data.h"
#include "linalg.h"
#include <iostream>
using namespace std;
/*
int  main (int  argc, char  *argv[])
{

    int numberAssets=83;
    int numberReturns=700;
    Lattice returnMatrix(numberAssets, Vector(numberReturns)); // a matrix to store the return data

    //read the data from the file and store it into the return matrix
    string fileName="asset_returns.csv";
    readData(returnMatrix,fileName);          // returnMatrix[i][j] stores the asset i, return j value

*//*    //example on how to calculate the average return
    double mean=0;
    for(int i=0;i<numberAssets;i++){
        mean=0;
        for(int j=0;j<numberReturns;j++)
        {
            double temp=returnMatrix[i][j];
//        cout << "Asset " << i << ", Return "<<i<<"="<< temp<<"\n";
            mean=mean+temp/numberReturns;
        }
    }
*//*
    return 0;
}*/


// Example usage
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
//    (A - I).prn();

    return 0;
}
