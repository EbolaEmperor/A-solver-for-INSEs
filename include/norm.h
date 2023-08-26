#ifndef _NORM_H_
#define _NORM_H_

#include <vector>
#include <cmath>
#include "matrix.h"

class Norm{
public:
    virtual double operator () (const ColVector &vec) const = 0;
};

class Norm_p : public Norm{
private:
    double p;
public:
    Norm_p(const int &p): p(p) {}
    double operator () (const ColVector &vec) const;
};

class Norm_inf : public Norm{
public:
    double operator () (const ColVector &vec) const;
};

#endif