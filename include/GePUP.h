#ifndef _GEPUP_H_
#define _GEPUP_H_

#include "amgSolver.h"
#include "sparseMatrix.h"
#include "matrix.h"
#include "norm.h"
#include "idpair.h"
#include "RKTable.h"
#include "TimeFunction2D.h"
#include "Field.h"
#include <vector>

class GePUP{
private:
    int M;
    double tEnd;
    double dH, dT;
    double nu;
    double eps = 1e-9;
    TimeFunction2D *g[2], *initial[2];
    // g: forcing term;   initial: initial condition;
    bool noForcingTerm;
    Field u, w;
    amgSolver poissonDirichlet, poissonNeumann, IRK;

    inline int idx(const int &i, const int &j) const;
    inline int idx(const idpair &i) const;
    inline bool inRange(const idpair &i) const;
    inline bool isGhost1(const idpair &i) const;
    inline bool isGhost2(const idpair &i) const;
    double value(const ColVector &phi, const idpair &i) const;

    ColVector L(const ColVector &phi) const;
    Field L(const Field &u) const;
    ColVector Gd(const ColVector &phi, const int &d) const;
    ColVector D(const Field &u) const;
    Field Proj(const Field &u) const;
    ColVector D(const TimeFunction2D *g, const double &t) const;
    Field Duu(const Field &u) const;

public:
    GePUP();
};

class GePUP_IMEX : public GePUP{
};

class GePUP_ERK : public GePUP{
};

#endif