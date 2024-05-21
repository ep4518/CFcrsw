#include "linalg.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
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
    int k = 0;
    double a, B;
    Matrix x(Lattice(rows, Vector(1, 0.5)));  // Initial guess (zero.5 vector)
    Matrix s = b - (*this) * x; // Initial residual
    Matrix sOld = s;
    Matrix p = s;  // Residual vector

    do {
        a = s.dot(s) / p.dot((*this) * p);
        x = x + (a * p);
        sOld = s;
        s = s - a * ((*this) * p);
        B = s.dot(s) / sOld.dot(sOld);
        p = s + B * p;
        k++;
    } while (s.dot(s) > 0.00001);

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
            result(i, j) = A(i, j) - B(i, j);
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

// vertical stack three matrices
Matrix vstack(const Matrix &A, const Matrix &B, const Matrix &C) {
    if (A.getColumns() != B.getColumns() || A.getColumns() != C.getColumns()) {
        throw std::invalid_argument("vstack requires matrices of the same width (h * w).");
    }
    int h = A.getRows() + B.getRows() + C.getRows();
    Matrix result(h, A.getColumns());

    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = A(i, j);
        }
    }
    for (int i = 0; i < B.getRows(); ++i) {
        for (int j = 0; j < B.getColumns(); ++j) {
            result(i + A.getRows(), j) = B(i, j);
        }
    }
    for (int i = 0; i < C.getRows(); ++i) {
        for (int j = 0; j < C.getColumns(); ++j) {
            result(i + A.getRows() + B.getRows(), j) = C(i, j);
        }
    }
    return result;
}

// horizontal stack three matrices
Matrix hstack(const Matrix &A, const Matrix &B, const Matrix &C) {
    if (A.getRows() != B.getRows() || A.getRows() != C.getRows()) {
        throw std::invalid_argument("hstack requires matrices of the same height (h * w).");
    }
    int w = A.getColumns() + B.getColumns() + C.getColumns();
    Matrix result(A.getRows(), w);
    // Copy the elements of A into the result matrix

    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = A(i, j);
        }
    }

    // Copy the elements of B into the result matrix
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < B.getColumns(); ++j) {
            result(i, j + A.getColumns()) = B(i, j);
        }
    }

    // Copy the elements of C into the result matrix
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < C.getColumns(); ++j) {
            result(i, j + A.getColumns() + B.getColumns()) = C(i, j);
        }
    }

    return result;
}