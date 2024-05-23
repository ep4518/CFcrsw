// main.cpp
#include "read_data.h"
#include "write_data.h"
#include "linalg.h"
#include "Markowitz.h"
#include <iostream>

using namespace std;

Matrix back_testing(const Matrix &optimal_weights, const Matrix &OOS_returns, const Matrix &target_returns);

int  main (int  argc, char  *argv[])
{

    int numberAssets=83;
    int numberReturns=700;
    Lattice returnMatrix(numberAssets, Vector(numberReturns)); // a matrix to store the return data

    //read the data from the file and store it into the return matrix
    string fileName="asset_returns.csv";
    readData(returnMatrix,fileName);          // returnMatrix[i][j] stores the asset i, return j value

    Matrix daily_returns(returnMatrix); // Assets - r * Days - c == 83 * 700

    Lattice tr =  {
            {0.   , 0.005, 0.01 , 0.015, 0.02 , 0.025, 0.03 , 0.035, 0.04 ,
            0.045, 0.05 , 0.055, 0.06 , 0.065, 0.07 , 0.075, 0.08 , 0.085,
            0.09 , 0.095, 0.1}
    };

    Matrix target_returns(tr);

    Result results[50];
    for (int i = 0; i < numberReturns - 100; i += 12) {
        int index = int(i/12);
        cout << "Moving to index " << index << endl;
        int start = index, mid = index + 100, end = index + 112;
        Matrix daily_returns_IS = daily_returns(0, numberAssets, start, mid);
        Matrix daily_returns_OOS = daily_returns(0, numberAssets, mid, end);
        Markowitz portfolio(daily_returns_IS, target_returns);
        portfolio.weights();
        Matrix df_optimal_weights = portfolio.weights();
        Matrix df_act_returns = back_testing(df_optimal_weights, daily_returns_OOS, target_returns);
        results[index] = {
                index,
                df_act_returns,
                df_optimal_weights
        };
    }

    results[0].back_test.prn();

    write_data(results, 50);

    return 0;
}

Matrix back_testing(const Matrix &optimal_weights, const Matrix &OOS_returns, const Matrix &target_returns) {
    int num_targ_rets = optimal_weights.getRows();
    int num_assets = OOS_returns.getRows();
    Markowitz OOS_rets(OOS_returns, target_returns);
    // Result matrix to store performance for each set of weights
    Matrix results(num_targ_rets, 3);

    for (int i = 0; i < num_targ_rets; i++) {
        // Extract the i-th row of optimal weights
        Matrix optimal_weights_row(num_assets, 1);
        for (int j = 0; j < num_assets; j++) {
            optimal_weights_row.insert(j, 0, optimal_weights(i, j));
        }
        Matrix act_ave_return = OOS_rets.mean().transpose() * optimal_weights_row;
        Matrix pf_cov = optimal_weights_row.transpose() * OOS_rets.cov() * optimal_weights_row;
        results.insert(i, 0, target_returns(0,i));
        results.insert(i, 1, act_ave_return(0,0));
        results.insert(i, 2, pf_cov(0,0));
    }
    return results;
}