#include "linalg.h"
#include <cmath>
#include <stdexcept>
//#include <cstddef>

Matrix Matrix::transpose() {
    Lattice transposed(columns, Vector(rows));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < columns; ++j) {
            transposed[j][i] = M[i][j];
        }
    }
    return {transposed};    // Clang-Tidy: Avoid repeating the return type from the declaration; use a braced initializer list instead
}

// Conjugate Gradient Solver
Matrix Matrix::solver(const Matrix &b) const {
    // requires non-singular matrix
    if (M.empty() || M[0].size() != M.size()) {
        throw std::invalid_argument("Matrix dimensions are not compatible for inversion");
    }
    // check b is a column vector
    if (b.getColumns() != 1 || b.getRows() != rows) {
        throw std::invalid_argument("b must be a column vector with the same number of rows as the matrix.");
    }
    // requires positive semi-definite matrix.
    // computationally expensive to calculate eigenvaules so checked in calculation    Matrix x(rows, 1);  // Initial guess (zero vector)

    Matrix x(rows, 1);  // Initial guess (zero vector)
    Matrix r = b;  // Residual vector
    Matrix p = r;  // Search direction
    Matrix Ap(rows, 1);
    double rsold = r.dot(r);

    for (int i = 0; i < rows; ++i) {
        Ap = (*this) * p;
        double alpha = rsold / p.dot(Ap);
        x = x + (alpha * p);
        r = r - (alpha * Ap);
        double rsnew = r.dot(r);
        if (std::sqrt(rsnew) < 1e-10) {
            break;
        }
        p = r + ((rsnew / rsold) * p);
        rsold = rsnew;
    }

    return x;  // Return result as a column Matrix
}

// Matrix multiplication
Matrix operator*(const Matrix& A, const Matrix& B) {
    if (A.getColumns() != B.getRows()) {
        throw std::invalid_argument("Matrix dimensions must agree for multiplication.");
    }
    Matrix result(A.getRows(), B.getColumns());
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < B.getColumns(); ++j) {
            for (int k = 0; k < A.getColumns(); ++k) {
                result(i, j) += A(i, k) * B(k, j);
            }
        }
    }
    return result;
}

// Matrix addition/ subtraction
Matrix operator+(const Matrix& A, const Matrix& B) {
    if (A.getRows() != B.getRows() || A.getColumns() != B.getColumns()) {
        throw std::invalid_argument("Matrix dimensions must agree for addition.");
    }
    Matrix result(A.getRows(), A.getColumns());
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = A(i, j) + B(i, j);
        }
    }
    return result;
}

Matrix operator-(const Matrix& A, const Matrix& B) {
    if (A.getRows() != B.getRows() || A.getColumns() != B.getColumns()) {
        throw std::invalid_argument("Matrix dimensions must agree for addition.");
    }
    Matrix result(A.getRows(), A.getColumns());
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = A(i, j) + B(i, j);
        }
    }
    return result;
}

// Scalar multiplication
Matrix operator*(const double& a, const Matrix& A) {
    Matrix result(A.getRows(), A.getColumns());
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = a * A(i, j);
        }
    }
    return result;
}

// Dot product of two column matrices
double Matrix::dot(const Matrix& A) const {
    if (columns != 1 || A.getColumns() != 1 || rows != A.getRows()) {
        throw std::invalid_argument("Dot product requires column matrices of the same size.");
    }
    double result = 0.0;
    for (int i = 0; i < rows; ++i) {
        result += M[i][0] * A(i, 0);
    }
    return result;
}