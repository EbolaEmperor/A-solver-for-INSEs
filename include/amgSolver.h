#ifndef _AMG_SOLVER_H_
#define _AMG_SOLVER_H_

#include "sparseMatrix.h"
#include "matrix.h"
#include <vector>
#include <cstring>

class amgSolver{
private:
    int AMG_DIRECT_N = 128;
    int AMG_SMOOTH_ITER = 3;
    double AMG_STRONG_THERESHOLD = 0.2;
    bool pureNeumann = false;
    std::vector<SparseMatrix> Ah, Ph, Rh;
    // Ah:算子  Ph:插值  Rh:限制

    ColVector VC(const int &d, ColVector x, const ColVector &rhs) const;
    ColVector FMG(const int &d, const ColVector &rhs) const;
    void generateGrid(const int &d);
    SparseMatrix getInterpolator(const SparseMatrix &A);
    std::vector<int> getCorsetPoints(const SparseMatrix &A);

public:
    ~amgSolver();
    void setPureNeumann();
    void setSmoothIter(const int &n);
    void setDirectN(const int &n);
    void setStrongThereshold(const double &d);

    void generateGrid(const SparseMatrix &A);
    ColVector solve(const ColVector &rhs, const std::string &method, const int &maxIter, const double &eps) const;
};

#endif