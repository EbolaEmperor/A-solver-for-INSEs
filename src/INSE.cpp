#include "INSE.h"
#include "RKTable.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

//------------------------------------------------------Settings----------------------------------------------------------------

INSE_Solver::INSE_Solver(){
    M = 0;
    dH = dT = 0.0;
    g[0] = g[1] = nullptr;
    initial[0] = initial[1] = nullptr;
    noForcingTerm = false;
}


void INSE_Solver::setNoForcingTerm(){
    noForcingTerm = true;
}

void INSE_Solver::setForcingTerm(TimeFunction2D *_g[2]){
    g[0] = _g[0];
    g[1] = _g[1];
}

void INSE_Solver::setInitial(TimeFunction2D *_initial[2]){
    initial[0] = _initial[0];
    initial[1] = _initial[1];
}

void INSE_Solver::setGridSize(const int &_M){
    M = _M;
    dH = 1.0 / M;
}

void INSE_Solver::setEndTime(const double &_tEnd){
    tEnd = _tEnd;
}

void INSE_Solver::setReynolds(const double &_R){
    nu = 1.0/_R;
}

void INSE_Solver::setNu(const double &_nu){
    nu = _nu;
}

void INSE_Solver::setEps(const double &_eps){
    eps = _eps;
}

void INSE_Solver::setTimeStepWithCaurant(const double &caurant, const double &maxux, const double &maxuy){
    if(dH==0.0){
        std::cerr << "[Error] setTimeStepWithCaurant: Please set grid size first." << std::endl;
    } else {
        dT = caurant / (maxux/dH + maxuy/dH);
    }
}

//---------------------------------------------------Discrete Operators---------------------------------------------------------

ColVector INSE_Solver::L(const ColVector &phi){
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            double t = 0.0;
            t += -60 * phi(idx(i,j));
            t += 16 * (phi(idx(i+1,j)) + phi(idx(i-1,j)) + phi(idx(i,j+1)) + phi(idx(i,j-1)));
            t += -(phi(idx(i+2,j)) + phi(idx(i-2,j)) + phi(idx(i,j+2)) + phi(idx(i,j-2)));
            res(idx(i,j)) = t / (12.0*dH*dH);
        }
    return res;
}

ColVector INSE_Solver::Duu(const std::vector<ColVector> &u, const int &k){
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) = ( F(u[0],u[k], idpair(i,j), idpair(1,0)) - F(u[0],u[k], idpair(i,j), idpair(-1,0))
                       + F(u[1],u[k], idpair(i,j), idpair(0,1)) - F(u[1],u[k], idpair(i,j), idpair(0,-1)) ) / dH;
        }
    return res;
}

ColVector INSE_Solver::XE(const std::vector<ColVector> &u, const int &k, const double &t){
    ColVector res = -Duu(u,k);
    if(noForcingTerm) return res;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) += g[k]->int2D_order6(i*dH, (i+1)*dH, j*dH, (j+1)*dH, t) * M * M;
        }
    return res;
}

double INSE_Solver::face(const ColVector &phi, const idpair &i, const idpair &ed){
    return 7.0/12.0 * (phi(idx(i)) + phi(idx(i+ed))) - 1.0/12.0 * (phi(idx(i-ed)) + phi(idx(i+2*ed)));
}

double INSE_Solver::GdVertical(const ColVector &phi, const idpair &i, const idpair &ed){
    idpair edp = ed.transDirection();
    return (face(phi,i+edp,ed) - face(phi,i-edp,ed)) / (2*dH);
}

double INSE_Solver::F(const ColVector &phi, const ColVector &psi, const idpair &i, const idpair &ed){
    return face(phi,i,ed)*face(psi,i,ed) + dH*dH/12.0 * GdVertical(phi,i,ed)*GdVertical(psi,i,ed);
}

ColVector INSE_Solver::D(const std::vector<ColVector> &u){
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) = ( -u[0](idx(i+2,j)) + 8.0*u[0](idx(i+1,j)) - 8.0*u[0](idx(i-1,j)) + u[0](idx(i-2,j))
                              -u[1](idx(i,j+2)) + 8.0*u[1](idx(i,j+1)) - 8.0*u[1](idx(i,j-1)) + u[1](idx(i,j-2)) ) / (12.0*dH);
        }
    return res;
}

ColVector INSE_Solver::Gd(const ColVector &phi, const int &d){
    ColVector res(M*M);
    idpair ed;
    ed[d] = 1;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            idpair ib(i,j);
            res(idx(ib)) = ( -phi(idx(ib+2*ed)) + 8.0*phi(idx(ib+ed)) - 8.0*phi(idx(ib-ed)) + phi(idx(ib-2*ed)) ) / (12.0*dH);
        }
    return res;
}

std::vector<ColVector> INSE_Solver::Proj(const std::vector<ColVector> &u){
    ColVector rhs = D(u);
    ColVector tmp = invLsolver.solve(rhs, "FMG", 7, eps);
    std::vector<ColVector> res;
    res.push_back(u[0] + Gd(tmp,0));
    res.push_back(u[1] + Gd(tmp,1));
    return res;
}

//---------------------------------------------------Sol Value and Index--------------------------------------------------------

int INSE_Solver::idx(const int &i, const int &j){
    return (i+M)%M * M + (j+M)%M;
}

int INSE_Solver::idx(const idpair &x){
    return (x[0]+M)%M * M + (x[1]+M)%M;
}

//---------------------------------------------------output and check error---------------------------------------------------------

void INSE_Solver::output(const std::string &outname){
    std::cout << "--------------------------------------------------------------" << std::endl;
    std::ofstream out(outname);
    out << std::fixed << std::setprecision(16);
    out << sol_u[0].T() << std::endl;
    out << sol_u[1].T() << std::endl;
    out << sol_p.T() << std::endl;
    out.close();
    std::cout << "Result has been saved to " << outname << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
}

std::vector<double> INSE_Solver::checkerr(TimeFunction2D *truesol_u[2], TimeFunction2D *truesol_p, const Norm &norm){
    ColVector trueres[3];
    for(int i = 0; i < 3; i++)
        trueres[i] = ColVector(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            trueres[0](idx(i,j)) = truesol_u[0]->accInt2D(i*dH, (i+1)*dH, j*dH, (j+1)*dH, tEnd) / (dH*dH);
            trueres[1](idx(i,j)) = truesol_u[1]->accInt2D(i*dH, (i+1)*dH, j*dH, (j+1)*dH, tEnd) / (dH*dH);
            trueres[2](idx(i,j)) = truesol_p->accInt2D(i*dH, (i+1)*dH, j*dH, (j+1)*dH, tEnd) / (dH*dH);
        }
    std::vector<double> res;
    res.push_back(norm(trueres[0]-sol_u[0]));
    res.push_back(norm(trueres[1]-sol_u[1]));
    res.push_back(norm(trueres[2]-sol_p));
    return res;
}

//-----------------------------------------------------Solve Process----------------------------------------------------------

void INSE_Solver::solve(){
    using namespace ERK_ESDIRK_Table;

    std::cout << "Setting initial values..." << std::endl;
    for(int d = 0; d < 2; d++){
        sol_u.push_back(ColVector(M*M));
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++){
                sol_u[d](idx(i,j)) = initial[d]->accInt2D(i*dH, (i+1)*dH, j*dH, (j+1)*dH, 0) * M * M;
            }
    }
    //std::cout << sol_u[0].T() << std::endl;
    
    std::vector<ColVector> tmpsol[stage], XEtmpsol[stage], PXEtmpsol[stage], Ltmpsol[stage];
    for(int s = 0; s < stage; s++){
        for(int d = 0; d < 2; d++){
            tmpsol[s].push_back(ColVector());
            XEtmpsol[s].push_back(ColVector());
            PXEtmpsol[s].push_back(ColVector());
            Ltmpsol[s].push_back(ColVector());
        }
    }

    // Initialize the solver for inverse L operator
    matrixElements.clear();
    double con1 = 5.0/(dH*dH);
    double con2 = -4.0/(3.0*dH*dH);
    double con3 = 1.0/(12.0*dH*dH);
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            int r = i * M + j;
            matrixElements.push_back( Triple(r, idx(i,j), con1) );
            matrixElements.push_back( Triple(r, idx(i+1,j), con2) );
            matrixElements.push_back( Triple(r, idx(i-1,j), con2) );
            matrixElements.push_back( Triple(r, idx(i,j+1), con2) );
            matrixElements.push_back( Triple(r, idx(i,j-1), con2) );
            matrixElements.push_back( Triple(r, idx(i+2,j), con3) );
            matrixElements.push_back( Triple(r, idx(i-2,j), con3) );
            matrixElements.push_back( Triple(r, idx(i,j+2), con3) );
            matrixElements.push_back( Triple(r, idx(i,j-2), con3) );
        }
    }
    invLsolver.generateGrid(SparseMatrix(M*M, M*M, matrixElements));
    invLsolver.setPureNeumann();
    matrixElements.clear();

    // Initialize the solver for Implicit RK
    con1 = 1.0 + dT*nu*gma*5.0/(dH*dH);
    con2 = -dT*nu*gma*4.0/(3.0*dH*dH);
    con3 = dT*nu*gma/(12.0*dH*dH);
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            int r = i * M + j;
            matrixElements.push_back( Triple(r, idx(i,j), con1) );
            matrixElements.push_back( Triple(r, idx(i+1,j), con2) );
            matrixElements.push_back( Triple(r, idx(i-1,j), con2) );
            matrixElements.push_back( Triple(r, idx(i,j+1), con2) );
            matrixElements.push_back( Triple(r, idx(i,j-1), con2) );
            matrixElements.push_back( Triple(r, idx(i+2,j), con3) );
            matrixElements.push_back( Triple(r, idx(i-2,j), con3) );
            matrixElements.push_back( Triple(r, idx(i,j+2), con3) );
            matrixElements.push_back( Triple(r, idx(i,j-2), con3) );
        }
    }
    solver.setStrongThereshold(std::min(0.2, 5*nu));
    solver.generateGrid(SparseMatrix(M*M, M*M, matrixElements));
    matrixElements.clear();

    for(double t = 0.0; t+1e-12 < tEnd; t += dT){
        // multi stage
        auto cl = clock();
        std::cout << "Time: " << t << ". ";
        for(int d = 0; d < 2; d++){
            tmpsol[0][d] = sol_u[d];
            XEtmpsol[0][d] = XE(sol_u, d, t);
            Ltmpsol[0][d] = L(sol_u[d]);
        }
        PXEtmpsol[0] = Proj(XEtmpsol[0]);

        for(int s = 1; s < stage; s++){
            double ts = t + c[s]*dT;
            for(int d = 0; d < 2; d++){
                ColVector rhs = sol_u[d];
                for(int j = 0; j < s; j++){
                    rhs = rhs + dT*aE[s][j]*PXEtmpsol[j][d] + dT*nu*aI[s][j]*Ltmpsol[j][d];
                }
                tmpsol[s][d] = solver.solve(rhs, "FMG", 8, eps);
            }
            for(int d = 0; d < 2; d++){
                XEtmpsol[s][d] = XE(tmpsol[s], d, ts);
                Ltmpsol[s][d] = L(tmpsol[s][d]);
            }
            if(s < stage-1)
                PXEtmpsol[s] = Proj(XEtmpsol[s]);
        }

        // combine the results
        for(int d = 0; d < 2; d++){
            for(int s = 0; s < stage; s++){
                sol_u[d] = sol_u[d] + dT * b[s] * (XEtmpsol[s][d] + nu * Ltmpsol[s][d]);
            }
        }
        sol_u = Proj(sol_u);
        //std::cout << sol_u[0].T() << std::endl;
        std::cout << "Solved in " << (double)(clock()-cl)/CLOCKS_PER_SEC << "s" << std::endl;
    }
    
    std::cout << "The velocity field is solved. Start solving the pressure..." << std::endl;
    std::vector<ColVector> tmp;
    for(int d = 0; d < 2; d++){
        tmp.push_back(-Duu(sol_u, d) + nu * L(sol_u[d]));
    }
    if(!noForcingTerm){
        for(int d = 0; d < 2; d++)
            for(int i = 0; i < M; i++)
                for(int j = 0; j < M; j++){
                    tmp[d](idx(i,j)) += g[d]->accInt2D_order6(i*dH, (i+1)*dH, j*dH, (j+1)*dH, 0) * M * M;
                }
    }
    sol_p = invLsolver.solve(-D(tmp), "FMG", 7, eps);

    std::cout << "All solved." << std::endl;
}