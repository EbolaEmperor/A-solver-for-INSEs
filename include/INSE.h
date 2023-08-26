#ifndef _INSE_H_
#define _INSE_H_

#include "matrix.h"
#include "sparseMatrix.h"
#include "amgSolver.h"
#include "norm.h"
#include "TimeFunction2D.h"
#include "idpair.h"
#include <cstring>

class INSE_Solver{
private:
    int M;
    double tEnd;
    double dH, dT;
    double nu;
    double eps = 1e-9;
    amgSolver invLsolver, solver;
    TimeFunction2D *g[2], *initial[2];
    // g: forcing term;   initial: initial condition;
    bool noForcingTerm;
    ColVector sol_p;
    std::vector<ColVector> sol_u;
    std::vector<Triple> matrixElements;

    // Discrete operators
    ColVector L(const ColVector &phi);
    ColVector XE(const std::vector<ColVector> &u, const int &d, const double &t);
    ColVector Duu(const std::vector<ColVector> &u, const int &d);
    double F(const ColVector &phi, const ColVector &psi, const idpair &i, const idpair &ed);
    double face(const ColVector &phi, const idpair &i, const idpair &ed);
    double GdVertical(const ColVector &phi, const idpair &i, const idpair &ed);
    ColVector D(const std::vector<ColVector> &u);
    ColVector Gd(const ColVector &phi, const int &d);
    std::vector<ColVector> Proj(const std::vector<ColVector> &u);

    // if i and j are both in range, return the corresponding value,
    // else return the ghost cell value with corresponding bondaries.
    double solValue(const ColVector &phi, const int &i, const int &j, const double &t);
    void setMatrix(const int &r, const int &i, const int &j, const double &k);

    // map the 2D index into 1D index
    int idx(const int &i, const int &j);
    int idx(const idpair &x);

public:
    INSE_Solver();
    void solve();
    void setNoForcingTerm();
    void setForcingTerm(TimeFunction2D *_g[2]);
    void setInitial(TimeFunction2D *_initial[2]);
    void setGridSize(const int &_M);
    void setEndTime(const double &_tEnd);
    void setReynolds(const double &_R);
    void setNu(const double &_nu);
    void setEps(const double &_eps);
    void setTimeStepWithCaurant(const double &caurant, const double &maxux, const double &maxuy);
    void output(const std::string &outname);
    std::vector<double> checkerr(TimeFunction2D *truesol_u[2], TimeFunction2D *truesol_p, const Norm &norm);
};

#endif