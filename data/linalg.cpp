#include "linalg.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

void Matrix::prn() const {
    std::cout << "[\n";
    for (const auto& row : M) {
        std::cout << " [ ";
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "]\n";
}

Matrix Matrix::operator()(int row_start, int row_end, int col_start, int col_end) const {
    if (row_start < 0 || row_end > rows || col_start < 0 || col_end > columns) {
        throw std::out_of_range("Index out of bounds");
    }
    Lattice result;
    for (int i = row_start; i < row_end; ++i) {
        Vector row;
        for (int j = col_start; j < col_end; ++j) {
            row.push_back(M[i][j]);
        }
        result.push_back(row);
    }
    return Matrix(result);
}

// Transpose of the matrix
Matrix Matrix::transpose() const {
    Lattice transposed(columns, Vector(rows));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < columns; ++j) {
            transposed[j][i] = M[i][j];
        }
    }
    return {transposed};
}

// Conjugate Gradient Solver
Matrix Matrix::solver(const Matrix &b) const {
    if (M.empty() || rows != columns) {
        throw std::invalid_argument("Matrix dimensions are not compatible for inversion");
    }
    if (b.getColumns() != 1 || b.getRows() != rows) {
        throw std::invalid_argument("b must be a column vector with the same number of rows as the matrix.");
    }

    Matrix x(Lattice(rows, Vector(1, 1.0)));
    Matrix r = b - (*this) * x;
    Matrix p = r;
    double rsold = r.dot(r);
    double rsnew;

    for (int k = 0; k < b.getRows(); ++k) {
        Matrix Ap = (*this) * p;
        double alpha = rsold / p.dot(Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r.dot(r);
        if (std::sqrt(rsnew) < 1e-7) break;
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    }

    return x;
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

// Matrix addition
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

// Matrix subtraction
Matrix operator-(const Matrix& A, const Matrix& B) {
    if (A.getRows() != B.getRows() || A.getColumns() != B.getColumns()) {
        throw std::invalid_argument("Matrix dimensions must agree for subtraction.");
    }
    Matrix result(A.getRows(), A.getColumns());
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = A(i, j) - B(i, j);
        }
    }
    return result;
}

// Unary - operator for matrix
Matrix Matrix::operator-() const {
    Matrix result(rows, columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            result(i, j) = -M[i][j];
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

// LU Decomposition Solver
void Matrix::luDecompose(Matrix& L, Matrix& U) const {
    if (rows != columns) {
        throw std::invalid_argument("Matrix must be square for LU decomposition.");
    }
    int n = rows;
    L = Matrix(n, n);
    U = *this;

    for (int i = 0; i < n; ++i) {
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += (L(i, j) * U(j, k));
            }
            U(i, k) -= sum;
        }
        for (int k = i; k < n; ++k) {
            if (i == k) {
                L(i, i) = 1.0;
            } else {
                double sum = 0.0;
                for (int j = 0; j < i; ++j) {
                    sum += (L(k, j) * U(j, i));
                }
                L(k, i) = (U(k, i) - sum) / U(i, i);
            }
        }
    }
}

Matrix Matrix::forwardSubstitution(const Matrix& L, const Matrix& b) const {
    int n = L.getRows();
    Matrix y(n, 1);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L(i, j) * y(j, 0);
        }
        y(i, 0) = (b(i, 0) - sum) / L(i, i);
    }
    return y;
}

Matrix Matrix::backwardSubstitution(const Matrix& U, const Matrix& y) const {
    int n = U.getRows();
    Matrix x(n, 1);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += U(i, j) * x(j, 0);
        }
        x(i, 0) = (y(i, 0) - sum) / U(i, i);
    }
    return x;
}

Matrix Matrix::solveLU(const Matrix& b) const {
    if (b.getColumns() != 1 || b.getRows() != rows) {
        throw std::invalid_argument("b must be a column vector with the same number of rows as the matrix.");
    }

    Matrix L, U;
    luDecompose(L, U);

    Matrix y = forwardSubstitution(L, b);
    Matrix x = backwardSubstitution(U, y);
    return x;
}

// Vertical stack of three matrices
Matrix vstack(const Matrix &A, const Matrix &B, const Matrix &C) {
    if (A.getColumns() != B.getColumns() || A.getColumns() != C.getColumns()) {
        throw std::invalid_argument("vstack requires matrices of the same width.");
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

// Horizontal stack of three matrices
Matrix hstack(const Matrix &A, const Matrix &B, const Matrix &C) {
    if (A.getRows() != B.getRows() || A.getRows() != C.getRows()) {
        throw std::invalid_argument("hstack requires matrices of the same height.");
    }
    int w = A.getColumns() + B.getColumns() + C.getColumns();
    Matrix result(A.getRows(), w);

    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getColumns(); ++j) {
            result(i, j) = A(i, j);
        }
    }
    for (int i = 0; i < B.getRows(); ++i) {
        for (int j = 0; j < B.getColumns(); ++j) {
            result(i, j + A.getColumns()) = B(i, j);
        }
    }
    for (int i = 0; i < C.getRows(); ++i) {
        for (int j = 0; j < C.getColumns(); ++j) {
            result(i, j + A.getColumns() + B.getColumns()) = C(i, j);
        }
    }
    return result;
}
