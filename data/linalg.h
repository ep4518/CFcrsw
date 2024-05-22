// linalg.h
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
    // Default constructor
    Matrix() : rows(0), columns(0), M() {}

    // Construct from pre-existing Matrix
    Matrix(Lattice _M) {M = _M; rows = _M.size(); columns = _M[0].size();}

    // Construct an empty matrix with given dimensions
    Matrix(int _rows, int _columns) : rows(_rows), columns(_columns), M(_rows, Vector(_columns, 0.0)) {}

    // Access element (const version)
    const double& operator()(size_t i, size_t j) const { return M[i][j]; }

    // Method to insert an element into the matrix
    void insert(size_t i, size_t j, double value) {
        if (i >= rows || j >= columns) {
            throw std::out_of_range("Matrix indices out of range");
        }
        M[i][j] = value;
    }

    // pandas like shape method
    void shape() {printf("(%d, %d)\n",rows, columns);}

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

    // Unary - operator for matrix
    Matrix operator-() const {
        Matrix result(rows, columns);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                result(i, j) = -M[i][j];
            }
        }
        return result;
    }

    // Splicing operator
    Matrix operator()(int row_start, int row_end, int col_start, int col_end) const {
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

    Matrix transpose();

    // Skip inversion. Utilise conjugate gradient algorithm (see .pdf/ wiki)
    // Mx = b => x*
    //  https://en.wikipedia.org/wiki/Conjugate_gradient_method
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

// hstack and vstack for constructing Q
Matrix vstack(const Matrix &A, const Matrix &B, const Matrix &C);
Matrix hstack(const Matrix &A, const Matrix &B, const Matrix &C);

#endif