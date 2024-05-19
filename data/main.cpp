#include "read_data.h"
#include "linalg.h"
#include <iostream>
using namespace std;

int  main (int  argc, char  *argv[])
{

    int numberAssets=83;
    int numberReturns=700;
    Matrix returnMatrix(numberAssets, Vector(numberReturns)); // a matrix to store the return data

    //read the data from the file and store it into the return matrix
    string fileName="asset_returns.csv";
    readData(returnMatrix,fileName);          // returnMatrix[i][j] stores the asset i, return j value





/*    //example on how to calculate the average return
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
*/
    return 0;
}