#include "TimeFunction2D.h"
#include <cmath>

double TimeFunction2D::intFixX(const double &x, const double &d, const double &u, const double &t) const{
    return (u-d)/6.0 * (at(x,d,t) + 4.0 * at(x,(d+u)/2.0,t) + at(x,u,t));
}

double TimeFunction2D::intFixY(const double &y, const double &d, const double &u, const double &t) const{
    return (u-d)/6.0 * (at(d,y,t) + 4.0 * at((d+u)/2.0,y,t) + at(u,y,t));
}

double TimeFunction2D::int2D(const double &l, const double &r, const double &d, const double &u, const double &t) const{
    static const double simpsonCoef[3][3] = {
        1, 4, 1,
        4, 16, 4,
        1, 4, 1
    };
    const double dx = (r-l) / 2;
    const double dy = (u-d) / 2;
    double res = 0;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++){
            res += simpsonCoef[i][j] * at(l+i*dx, d+j*dy, t);
        }
    return res * dx * dy / 9;
}

double TimeFunction2D::_accInt2D(const double &l, const double &r, const double &d, const double &u, const double &t, const double &A) const{
    const double midx = (l+r)/2;
    const double midy = (d+u)/2;
    const double A1 = int2D(l,midx,d,midy,t);
    const double A2 = int2D(midx,r,d,midy,t);
    const double A3 = int2D(l,midx,midy,u,t);
    const double A4 = int2D(midx,r,midy,u,t);
    return fabs(A1+A2+A3+A4-A)<1e-14 ? A1+A2+A3+A4 : 
        _accInt2D(l,midx,d,midy,t,A1)+_accInt2D(midx,r,d,midy,t,A2)+
        _accInt2D(l,midx,midy,u,t,A3)+_accInt2D(midx,r,midy,u,t,A4);

}

double TimeFunction2D::accInt2D(const double &l, const double &r, const double &d, const double &u, const double &t) const{
    return _accInt2D(l, r, d, u, t, int2D(l,r,d,u,t));
}

double TimeFunction2D::intFixX_order6(const double &x, const double &d, const double &u, const double &t) const{
    double delta = u-d;
    return delta/90.0 * (7.0*at(x,d,t) + 32.0*at(x,d+0.25*delta,t) + 12.0*at(x,d+0.5*delta,t) + 32.0*at(x,d+0.75*delta,t) + 7.0*at(x,u,t));
}

double TimeFunction2D::intFixY_order6(const double &y, const double &d, const double &u, const double &t) const{
    double delta = u-d;
    return delta/90.0 * (7.0*at(d,y,t) + 32.0*at(d+0.25*delta,y,t) + 12.0*at(d+0.5*delta,y,t) + 32.0*at(d+0.75*delta,y,t) + 7.0*at(u,y,t));
}

double TimeFunction2D::int2D_order6(const double &l, const double &r, const double &d, const double &u, const double &t) const{
    static const double simpsonCoef[5][5] = {
        49, 224, 84, 224, 49, 
        224, 1024, 384, 1024, 224, 
        84, 384, 144, 384, 84, 
        224, 1024, 384, 1024, 224, 
        49, 224, 84, 224, 49
    };
    const double dx = (r-l) / 4;
    const double dy = (u-d) / 4;
    double res = 0;
    for(int i = 0; i < 5; i++)
        for(int j = 0; j < 5; j++){
            res += simpsonCoef[i][j] * at(l+i*dx, d+j*dy, t);
        }
    return res * (r-l) * (u-d) / 8100;
}

double TimeFunction2D::_accInt2D_order6(const double &l, const double &r, const double &d, const double &u, const double &t, const double &A) const{
    const double midx = (l+r)/2;
    const double midy = (d+u)/2;
    const double A1 = int2D_order6(l,midx,d,midy,t);
    const double A2 = int2D_order6(midx,r,d,midy,t);
    const double A3 = int2D_order6(l,midx,midy,u,t);
    const double A4 = int2D_order6(midx,r,midy,u,t);
    return fabs(A1+A2+A3+A4-A)<1e-14 ? A1+A2+A3+A4 : 
        _accInt2D_order6(l,midx,d,midy,t,A1)+_accInt2D_order6(midx,r,d,midy,t,A2)+
        _accInt2D_order6(l,midx,midy,u,t,A3)+_accInt2D_order6(midx,r,midy,u,t,A4);

}

double TimeFunction2D::accInt2D_order6(const double &l, const double &r, const double &d, const double &u, const double &t) const{
    return _accInt2D_order6(l, r, d, u, t, int2D_order6(l,r,d,u,t));
}