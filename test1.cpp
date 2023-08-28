#include <bits/stdc++.h>
#include "GePUP_IMEX.h"
using namespace std;

const double pi = acos(-1);

class Funcux : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return sin(pi*x)*sin(pi*x)*sin(2*pi*y) + 2*pi*cos(2*pi*x)*sin(2*pi*y);
    }
};

class Funcuy : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return -sin(pi*y)*sin(pi*y)*sin(2*pi*x) + 2*pi*sin(2*pi*x)*cos(2*pi*y);
    }
};

int main(int argc, char* argv[]){
    TimeFunction2D *u[2];
    u[0] = new Funcux();
    u[1] = new Funcuy();

    GePUP_IMEX solver;
    solver.setGridSize(stoi(argv[1]));
    solver.setEndTime(0.5);
    solver.setReynolds(1e4);
    solver.setTimeStepWithCaurant(stod(argv[2]), 2*pi+1, 2*pi+1);
    solver.setNoForcingTerm();
    solver.setInitial(u);
    solver.setEps(1e-9);
    solver.solve();
    solver.output("result.txt");
    return 0;
}