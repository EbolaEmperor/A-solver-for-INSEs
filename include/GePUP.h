#ifndef _GEPUP_H_
#define _GEPUP_H_

#include "amgSolver.h"
#include "sparseMatrix.h"
#include "matrix.h"
#include "norm.h"
#include "idpair.h"
#include "TimeFunction2D.h"
#include "Field.h"
#include <vector>

class GePUP{
protected:
    int M;
    double tEnd;
    double dH, dT;
    double nu;
    double eps = 1e-9;
    TimeFunction2D *g[2], *initial[2];
    // g: forcing term;   initial: initial condition;
    bool noForcingTerm = false;
    Field u, w;
    amgSolver poissonDirichlet, poissonNeumann;

    inline int idx(const int &i, const int &j) const {return i*M + j;}
    inline int idx(const idpair &i) const {return i[0]*M + i[1];}
    inline bool inRange(const idpair &i) const {return i[0]>=0 && i[0]<M && i[1]>=0 && i[1]<M;}
    inline bool isGhost1(const idpair &i) const {return i[0]==-1 || i[0]==M || i[1]==-1 || i[1]==M;}
    inline bool isGhost2(const idpair &i) const {return i[0]==-2 || i[0]==M+1 || i[1]==-2 || i[1]==M+1;}
    double value(const ColVector &phi, const idpair &i) const;

    ColVector L(const ColVector &phi) const;
    Field L(const Field &u) const;
    ColVector Gd(const ColVector &phi, const int &d) const;
    ColVector D(const Field &u) const;
    double D(const Field &u, const idpair &i) const;
    Field Proj(const Field &u) const;
    ColVector D(TimeFunction2D *const *, const double &t) const;
    Field Duu(const Field &u) const;
    Field XE(const Field &u, const double &t) const;
    ColVector bodyAve(TimeFunction2D *g, const double &t) const;
    double normalFace(TimeFunction2D *const *g, const double &t, const idpair &j, const idpair &ed) const;
    double normalPartialDivU(const Field &u, const idpair &j, const idpair &ed) const;
    double normalPartialQ(const Field &u, TimeFunction2D *const *g, const double &t, const idpair &j, const idpair &ed) const;
    double deltaUn(const Field &u, const idpair &j, const idpair &ed) const;
    void addRHSElementForNeumann(ColVector &rhs, const int &row, const idpair &i, const double &t, const double &coef) const;
    ColVector solveQ(const Field &u, TimeFunction2D *const *, const double &t) const;
    double face(const ColVector &phi, const idpair &i, const idpair &ed) const;
    double GdVer(const ColVector &phi, const idpair &i, const idpair &ed) const;
    double F(const ColVector &phi, const ColVector &psi, const idpair &i, const idpair &ed) const;

    void addElementForHomoDirichlet(std::vector<Triple> &elements, const int &row, const idpair &i, const double &coef);
    void addElementForNeumann(std::vector<Triple> &elements, const int &row, const idpair &i, const double &coef);
    void initialize();

public:
    void setNoForcingTerm();
    void setForcingTerm(TimeFunction2D *_g[2]);
    void setInitial(TimeFunction2D *_initial[2]);
    void setGridSize(const int &_M);
    void setEndTime(const double &_tEnd);
    void setReynolds(const double &_R);
    void setNu(const double &_nu);
    void setEps(const double &_eps);
    void setTimeStepWithCaurant(const double &caurant, const double &maxux, const double &maxuy);
};

#endif