//Markowitz.h
#ifndef MARKOWITZ_H
#define MARKOWITZ_H

#include "linalg.h"

#define TOLERANCE 1e-10
typedef enum {NODEBUG, DEBUG};

class Markowitz {
private:
    Matrix returns;
    Matrix target_returns;
    Matrix optimal_weights;
    int n;

public:
    Markowitz(const Matrix &_returns, const Matrix &_target_returns) {
        if (_target_returns.getRows() != 1) {
            throw std::invalid_argument("target returns should be a column vector");
        }
        returns = _returns; target_returns = _target_returns, n = _returns.getRows();
        this->weights();
    }

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

    Matrix getWeights() {return this->optimal_weights;}

    void NormTest();

};

#endif