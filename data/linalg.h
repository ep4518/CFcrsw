#ifndef _LINALG_H
#define _LINALG_H

#include <vector>
#include <iostream>

typedef std::vector<double> Vector;
typedef std::vector<Vector> Lattice;

class Matrix {
private:
    int rows, columns;
    Lattice M;
public:
    // Construct from pre-existing Matrix
    Matrix(Lattice _M) {M = _M; rows = _M.size(); columns = _M[0].size();}

    // Construct an empty matrix with given dimensions
    Matrix(int _rows, int _columns) : rows(_rows), columns(_columns), M(_rows, Vector(_columns, 0.0)) {}

    // Access element (const version)
    const double& operator()(size_t i, size_t j) const { return M[i][j]; }

    // Get number of rows
    int getRows() const { return rows; }

    // Get number of columns
    int getColumns() const { return columns; }

    // Get the underlying matrix
    const Lattice& getMatrix() const { return M; }

    // Access element (non-const version)
    double& operator()(size_t i, size_t j) {
        if (i >= rows || j >= columns) throw std::out_of_range("Matrix indices out of range");
        return M[i][j];
    }

    Matrix transpose();

    // Skip inversion. Utilise conjugate gradient algorithm (see .pdf/ wiki)
    // Mx = b => x*
    // https://en.wikipedia.org/wiki/Conjugate_gradient_method
    Matrix solver(const Matrix &b) const;

    // Matrix multiplication
    friend Matrix operator*(const Matrix& A, const Matrix& B);
    // Scalar multiplication
    friend Matrix operator*(const double& a, const Matrix&A);

    // Matrix addition/ subtraction
    friend Matrix operator+(const Matrix& A, const Matrix& B);
    friend Matrix operator-(const Matrix& A, const Matrix& B);

    // Dot product of two column matrices
    double dot(const Matrix& A) const;

    // Print the matrix
    void prn() const {
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
};

#endif