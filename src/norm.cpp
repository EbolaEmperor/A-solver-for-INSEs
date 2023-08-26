#include "norm.h"

double Norm_p::operator () (const ColVector &vec) const{
    double sum = 0;
    for(int i = 0; i < vec.size(); i++)
        sum += pow(fabs(vec(i)), p);
    sum /= vec.size();
    return pow(sum, 1.0/p);
}

double Norm_inf::operator () (const ColVector &vec) const{
    double mx = 0;
    for(int i = 0; i < vec.size(); i++)
        mx = std::max(mx, fabs(vec(i)));
    return mx;
}