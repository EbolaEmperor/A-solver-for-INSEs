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
};

#endif