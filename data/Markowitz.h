#ifndef MARKOWITZ_H
#define MARKOWITZ_H

#include "linalg.h"

class Markowitz {
private:
    Matrix returns;
    Matrix target_returns;
    int n;

public:
    Markowitz(const Matrix &_returns, const Matrix &_target_returns) {
        if (_target_returns.getRows() != 1) {
            throw std::invalid_argument("target returns should be a column vector");
        }
        returns = _returns; target_returns = _target_returns, n = _returns.getRows();}

    // Mean returns
    Matrix mean();

    // Covariance matrix
    Matrix cov();

    // Q matrix
    Matrix Q();

    // b Vector
    Matrix b(const double &target_return);

    // weights dataframe
    Matrix weights();

};

#endif