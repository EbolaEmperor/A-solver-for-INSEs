#ifndef _GEPUP_IMEX_H_
#define _GEPUP_IMEX_H_

#include "GePUP.h"

class GePUP_IMEX : public GePUP{
private:
    amgSolver IRK;
public:
    void initialize();
    void solve();
};

#endif