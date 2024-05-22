#include "Markowitz.h"
#include <iostream>

Matrix Markowitz::mean() {
    int assets = returns.getRows();
    int days = returns.getColumns();
    Matrix result(assets, 1);
    double mean;
    for (int i = 0; i < assets; i++) {
        mean = 0;
        for (int j = 0; j < days; j++)
            mean += returns(i, j);
        mean /= days;
        result.insert(i, 0, mean);
    }
    return result;
}

Matrix Markowitz::cov() {
    int assets = returns.getRows();
    int days = returns.getColumns();
    Matrix rBar = this->mean();
    Matrix result(assets, assets);
    double sum;
    for (int i = 0; i < assets; i++) {
        for (int j = 0; j < assets; j++) {
            sum = 0;
            for (int k = 0; k < days; k++) {
                sum += (returns(i, k) - rBar(i, 0)) * (returns(j, k) - rBar(j, 0));
            }
            sum /= (days - 1);
            result.insert(i, j, sum);
        }
    }
    return result;
}

Matrix Markowitz::Q() {
    Matrix e(Lattice(returns.getRows(), Vector(1, 1.0)));
    Matrix r = this->mean();
    Matrix zero(Lattice(1, Vector(1, 0.0)));
    Matrix top_row = hstack(this->cov(), -r, -e);
    Matrix middle_row = hstack(-r.transpose(), zero, zero);
    Matrix bottom_row = hstack(-e.transpose(), zero, zero);
    Matrix A = vstack(top_row, middle_row, bottom_row);
    return A;
}

Matrix Markowitz::b(const double &target_return) {
    Matrix zeros(Lattice(returns.getRows(), Vector(1, 0.0)));
    Matrix ret(Lattice(1, Vector(1, target_return)));
    Matrix one(Lattice(1, Vector(1, 1.0)));
    Matrix b = vstack(zeros, -ret, -one);
    return b;
}

Matrix Markowitz::weights() {
    int m = target_returns.getColumns();
    Matrix results(m, this->n + 2);

    for (int i = 0; i < m; i++) {
        Matrix x = this->Q().solver(this->b(target_returns(0, i)));
        for (int j = 0; j < this->n + 2; j++) {
            results.insert(i, j, x(j, 0));
        }
    }
    return results;
}

