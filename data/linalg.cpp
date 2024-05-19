#include "linalg.h"
#include <cmath>
#include <stdexcept>
//#include <cstddef>

Matrixx Matrixx::transpose() {
    Matrix transposed(columns, Vector(rows));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < columns; ++j) {
            transposed[j][i] = M[i][j];
        }
    }
    return {transposed};
}

/*Vector Matrixx::solver(Vector b) {
    // requires non-singular matrix
    if (M.empty() || M[0].size() != M.size()) {
        throw std::invalid_argument("Matrix dimensions are not compatible for inversion");
    }
    // requires positive semi-definite matrix.
    // computationally expensive to calculate eigenvaules so checked in calculation
    int k;
    double tol = 0.00001;  // Typical tolerance 10e − 6
    Vector x0, s0, p0;  // Needs to be defined <========================
    Vector x = x0, s = s0, p = p0, a, sOld, B;

    do {
        a = s.transpose() * s / p.transpose() * M * p;
        x = x + a * p;
        sOld = s;
        s = s - a * M * p;
        B = s.transpose() * s / sOld.transpose() * sOld;
        p = s + B * p;
        k+=1;
    } while (sOld.tranpsose()*sOld <= tol);    // Do we need a vector class that is child of matrix?

    return {x};
}*/

Matrixx operator*(const Matrix& A, const Matrix& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size()) {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
    }

    size_t rows = A.size();
    size_t cols = B[0].size();
    size_t common_dim = B.size();

    Matrix result(rows, std::vector<double>(cols, 0.0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < common_dim; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

Matrixx operator+(const Matrix& A, const Matrix& B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrix dimensions are not compatible for addition");
    }

    size_t rows = A.size();
    size_t cols = A[0].size();

    Matrix result(rows, std::vector<double>(cols));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

/*Vector operator*(const Matrix& C,const Vector& V)
{
    int d = C.size();
    Vector W(d);
    for (int j=0; j<d; j++)
    {
        W[j]=0.0;
        for (int l=0; l<d; l++) W[j]=W[j]+C[j][l]*V[l];
    }
    return W;
}

Vector operator+(const Vector& V,const Vector& W)
{
    int d = V.size();
    Vector U(d);
    for (int j=0; j<d; j++) U[j] = V[j] + W[j];
    return U;
}

Vector operator+(const double& a,const Vector& V)
{
    int d = V.size();
    Vector U(d);
    for (int j=0; j<d; j++) U[j] = a + V[j];
    return U;
}

Vector operator*(const double& a,const Vector& V)
{
    int d = V.size();
    Vector U(d);
    for (int j=0; j<d; j++) U[j] = a*V[j];
    return U;
}

Vector operator*(const Vector& V,const Vector& W)
{
    int d = V.size();
    Vector U(d);
    for (int j=0; j<d; j++) U[j] = V[j] * W[j];
    return U;
}

Vector exp(const Vector& V)
{
    int d = V.size();
    Vector U(d);
    for (int j=0; j<d; j++) U[j] = exp(V[j]);
    return U;
}

double operator^(const Vector& V,const Vector& W)
{
    double sum = 0.0;
    int d = V.size();
    for (int j=0; j<d; j++) sum = sum + V[j]*W[j];
    return sum;
}*/

