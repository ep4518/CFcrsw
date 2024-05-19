#ifndef _LINALG_H
#define _LINALG_H

#include <vector>

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;

class Matrixx {
private:
    int rows, columns;
    Matrix M;
public:
    // Construct from pre-existing Matrix
    Matrixx(Matrix _M) {M = _M; rows = _M.size(); columns = _M[0].size();}

    // Access element (const version)
    const double& operator()(size_t i, size_t j) const { return M[i][j]; }

    // Get number of rows
    int getRows() const { return rows; }

    // Get number of columns
    int getColumns() const { return columns; }

    // Get the underlying matrix
    const Matrix& getMatrix() const { return M; }

    Matrixx transpose();

    // Skip inversion. Utilise conjugate gradient algorithm (see .pdf/ wiki)
    // Mx = b => x*
    // https://en.wikipedia.org/wiki/Conjugate_gradient_method
    Vector solver(Vector);
};

//Vector operator*(const Matrix& C,const Vector& V);
//Vector operator*(const double& a,const Vector& V);
//Vector operator+(const double& a,const Vector& V);
//Vector operator+(const Vector& V,const Vector& W);
//Vector operator*(const Vector& V,const Vector& W);
//Vector exp(const Vector& V);
//double operator^(const Vector& V,const Vector& W);

Matrixx operator*(const Matrix& A, const Matrix& B);
Matrixx operator+(const Matrix& A, const Matrix& B);

#endif