#ifndef _FIELD_H_
#define _FIELD_H_

#include "matrix.h"

class Field{
private:
    ColVector u[2];
public:
    Field(){}
    Field(const int &M){
        u[0] = ColVector(M*M);
        u[1] = ColVector(M*M);
    }
    Field(const ColVector &ux, const ColVector &uy){
        u[0] = ux;
        u[1] = uy;
    }
    Field(ColVector &&ux, ColVector &&uy){
        u[0] = std::move(ux);
        u[1] = std::move(uy);
    }
    const ColVector& operator [] (const int &i) const{
        return u[i];
    }
    ColVector& operator [] (const int &i){
        return u[i];
    }
    Field& operator += (const Field &rhs){
        u[0] += rhs[0];
        u[1] += rhs[1];
        return (*this);
    }
    Field& operator -= (const Field &rhs){
        u[0] -= rhs[0];
        u[1] -= rhs[1];
        return (*this);
    }
    friend Field operator * (const double &k, const Field &rhs){
        Field res;
        res[0] = k*rhs[0];
        res[1] = k*rhs[1];
        return res;
    }
};

#endif