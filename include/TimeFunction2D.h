#ifndef _TIME_FUNCTION_2D_H_
#define _TIME_FUNCTION_2D_H_

class TimeFunction2D{
public:
    virtual double at (const double &x, const double &y, const double &t) const = 0;
    virtual double intFixX(const double &x, const double &d, const double &u, const double &t) const;
    virtual double intFixY(const double &y, const double &d, const double &u, const double &t) const;
    virtual double int2D(const double &l, const double &r, const double &d, const double &u, const double &t) const;
    virtual double _accInt2D(const double &l, const double &r, const double &d, const double &u, const double &t, const double &A) const;
    virtual double accInt2D(const double &l, const double &r, const double &d, const double &u, const double &t) const;
    
    virtual double intFixX_order6(const double &x, const double &d, const double &u, const double &t) const;
    virtual double intFixY_order6(const double &y, const double &d, const double &u, const double &t) const;
    virtual double int2D_order6(const double &l, const double &r, const double &d, const double &u, const double &t) const;
    virtual double _accInt2D_order6(const double &l, const double &r, const double &d, const double &u, const double &t, const double &A) const;
    virtual double accInt2D_order6(const double &l, const double &r, const double &d, const double &u, const double &t) const;
};

#endif