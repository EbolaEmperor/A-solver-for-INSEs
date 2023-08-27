#include "GePUP.h"

int GePUP::idx(const int &i, const int &j) const{
    return i*M + j;
}

int GePUP::idx(const idpair &i) const{
    return i[0]*M + i[1];
}

inline bool GePUP::inRange(const idpair &i) const{
    return i[0]>=0 && i[0]<M && i[1]>=0 && i[1]<M;
}

inline bool GePUP::isGhost1(const idpair &i) const{
    return i[0]==-1 || i[0]==M || i[1]==-1 || i[1]==M;
}

inline bool GePUP::isGhost2(const idpair &i) const{
    return i[0]==-2 || i[0]==M+1 || i[1]==-2 || i[1]==M+1;
}

double GePUP::value(const ColVector &phi, const idpair &i) const{
    if(inRange(i)) return phi(idx(i));
    else if(isGhost1(i)){
        idpair j=i, ed;
        if(i[0]==-1) j[0]=0, ed[0]=-1;
        if(i[0]==M) j[0]=M-1, ed[0]=1;
        if(i[1]==-1) j[1]=0, ed[1]=-1;
        if(i[1]==M) j[1]=M-1, ed[1]=1;
        return (-77.0*phi(idx(j)) + 43.0*phi(idx(j-ed)) - 17.0*phi(idx(j-2*ed)) + 3.0*phi(idx(j-3*ed))) / 12.0;
    } else if(isGhost2(i)){
        idpair j=i, ed;
        if(i[0]==-2) j[0]=0, ed[0]=-1;
        if(i[0]==M+1) j[0]=M-1, ed[0]=1;
        if(i[1]==-2) j[1]=0, ed[1]=-1;
        if(i[1]==M+1) j[1]=M-1, ed[1]=1;
        return (-505.0*phi(idx(j)) + 335.0*phi(idx(j-ed)) - 145.0*phi(idx(j-2*ed)) + 27.0*phi(idx(j-3*ed))) / 12.0;
    } else {
        std::cerr << "[Error] value:: out of range!" << std::endl;
        exit(-1);
        return -1;
    }
}

ColVector GePUP::L(const ColVector &phi) const{
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            idpair id(i,j);
            for(int d = 0; d < 2; d++){
                idpair ed; ed[d]=1;
                res(idx(id)) += -value(phi,id+2*ed) + 16.0*value(phi,id+ed) - 30.0*value(phi,id) 
                                + 16.0*value(phi,id-ed) - value(phi,id-2*ed);
            }
            res(idx(id)) /= 12.0*dH*dH;
        }
    return res;
}

Field GePUP::L(const Field &u) const{
    return Field(L(u[0]), L(u[1]));
}

ColVector GePUP::Gd(const ColVector &phi, const int &d) const{
    ColVector res(M*M);
    idpair ed; ed[d]=1;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            idpair id(i,j);
            res(idx(id)) = (-value(u[d],id+2*ed) + 8.0*value(u[d],id+ed) -8.0*value(u[d],id-ed) + value(u[d],id-2*ed)) / (12.0*dH);
        }
    return res;
}

ColVector GePUP::D(const Field &u) const{
    return Gd(u[0],0) + Gd(u[1],1);
}

Field GePUP::Proj(const Field &u) const{
    auto tmp = poissonDirichlet.solve(D(u), "FMG", 7, eps);
    return Field( u[0]-Gd(tmp,0), u[1]-Gd(tmp,1) );
}

ColVector GePUP::D(TimeFunction2D *const *g, const double &t) const{
    // Compute the body-average integration of div(g), with the divergence theorem.
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) = (- g[1]->intFixY_order6(j*dH, i*dH, (i+1)*dH, t)
                             + g[1]->intFixY_order6((j+1)*dH, i*dH, (i+1)*dH, t)
                             - g[0]->intFixX_order6(i*dH, j*dH, (j+1)*dH, t)
                             + g[0]->intFixX_order6((i+1)*dH, j*dH, (j+1)*dH, t)) / (dH*dH);
        }
    return res;
}

double GePUP::face(const ColVector &phi, const idpair &i, const idpair &ed) const{
    return (-value(phi,i+2*ed) + 7.0*value(phi,i+ed) + 7.0*value(phi,i) - value(phi,i-ed)) / 12.0;
}

double GePUP::GdVer(const ColVector &phi, const idpair &i, const idpair &ed) const{
    idpair edVer = ed.transDirection();
    return (face(phi,i+edVer,ed) - face(phi,i-edVer,ed)) / (2.0*dH);
}

double GePUP::F(const ColVector &phi, const ColVector &psi, const idpair &i, const idpair &ed) const{
    if(isGhost1(i+ed)) return 0;
    return face(phi,i,ed)*face(psi,i,ed) + (dH*dH/12.0) * GdVer(phi,i,ed)*GdVer(psi,i,ed);
}

Field GePUP::Duu(const Field &u) const{
    Field res(M);
    for(int dm = 0; dm < 2; dm++)
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++){
                idpair id(i,j);
                for(int d = 0; d < 2; d++){
                    idpair ed; ed[d]=1;
                    res[dm](idx(id)) += (F(u[d],u[dm],id,ed) - F(u[d],u[dm],id,ed)) / dH;
                }
            }
    return res;
}

ColVector GePUP::bodyAve(TimeFunction2D *g, const double &t) const{
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) = g->int2D_order6(i*dH, (i+1)*dH, j*dH, (j+1)*dH, t);
        }
    return res;
}

double GePUP::deltaUn(const Field &u, const idpair &j, const idpair &ed) const{
    int d = ed[0] ? 0 : 1;
    double fg = (ed[d]==1) ? 1 : -1;
    return fg/(48.0*dH*dH) * ( -755.0*u[d](idx(j)) + 493.0*u[d](idx(j-ed)) - 191.0*u[d](idx(j-2*ed)) + 33.0*u[d](idx(j-3*ed)) );
}

double GePUP::normalFace(TimeFunction2D *const *g, const double &t, const idpair &j, const idpair &ed) const{
}

double GePUP::normalPartialDivU(const Field &u, const idpair &j, const idpair &ed) const{
}

double GePUP::normalPartialQ(const Field &u, TimeFunction2D *const *g, const double &t, const idpair &j, const idpair &ed) const{
}

ColVector GePUP::solveQ(const Field &u, TimeFunction2D *const * g, const double &t) const{
    ColVector rhs = D(g,t) - D(Duu(u));
    // ...
    return poissonNeumann.solve(rhs, "FMG", 8, eps);
}

Field GePUP::XE(const Field &u, const double &t) const{
    Field res;
    res[0] = bodyAve(g[0],t);
    res[1] = bodyAve(g[1],t);
    res -= Duu(u);
    auto q = solveQ(u,g,t);
    res[0] -= Gd(q,0);
    res[1] -= Gd(q,1);
    return res;
}

void GePUP::addElementForHomoDirichlet(std::vector<Triple> &elements, const int &row, const idpair &i, const double &coef){
    if(inRange(i)) elements.emplace_back(row, idx(i), coef);
    else if(isGhost1(i)){
        idpair j=i, ed;
        if(i[0]==-1) j[0]=0, ed[0]=-1;
        if(i[0]==M) j[0]=M-1, ed[0]=1;
        if(i[1]==-1) j[1]=0, ed[1]=-1;
        if(i[1]==M) j[1]=M-1, ed[1]=1;
        elements.emplace_back(row, idx(j), -77.0/12.0*coef);
        elements.emplace_back(row, idx(j-ed), 43.0/12.0*coef);
        elements.emplace_back(row, idx(j-2*ed), -17.0/12.0*coef);
        elements.emplace_back(row, idx(j-3*ed), 3.0/12.0*coef);
    } else if(isGhost2(i)){
        idpair j=i, ed;
        if(i[0]==-2) j[0]=0, ed[0]=-1;
        if(i[0]==M+1) j[0]=M-1, ed[0]=1;
        if(i[1]==-2) j[1]=0, ed[1]=-1;
        if(i[1]==M+1) j[1]=M-1, ed[1]=1;
        elements.emplace_back(row, idx(j), -505.0/12.0*coef);
        elements.emplace_back(row, idx(j-ed), 335.0/12.0*coef);
        elements.emplace_back(row, idx(j-2*ed), -145.0/12.0*coef);
        elements.emplace_back(row, idx(j-3*ed), 27.0/12.0*coef);
    } else {
        std::cerr << "[Error] addElementForHomoDirichlet:: out of range!" << std::endl;
        exit(-1);
    }
}

void GePUP::addElementForNeumann(std::vector<Triple> &elements, const int &row, const idpair &i, const double &coef){
    if(inRange(i)) elements.emplace_back(row, idx(i), coef);
    else if(isGhost1(i)){
        idpair j=i, ed;
        if(i[0]==-1) j[0]=0, ed[0]=-1;
        if(i[0]==M) j[0]=M-1, ed[0]=1;
        if(i[1]==-1) j[1]=0, ed[1]=-1;
        if(i[1]==M) j[1]=M-1, ed[1]=1;
        elements.emplace_back(row, idx(j), 0.5*coef);
        elements.emplace_back(row, idx(j-ed), 0.9*coef);
        elements.emplace_back(row, idx(j-2*ed), -0.5*coef);
        elements.emplace_back(row, idx(j-3*ed), 0.1*coef);
    } else if(isGhost2(i)){
        idpair j=i, ed;
        if(i[0]==-2) j[0]=0, ed[0]=-1;
        if(i[0]==M+1) j[0]=M-1, ed[0]=1;
        if(i[1]==-2) j[1]=0, ed[1]=-1;
        if(i[1]==M+1) j[1]=M-1, ed[1]=1;
        elements.emplace_back(row, idx(j), -7.5*coef);
        elements.emplace_back(row, idx(j-ed), 14.5*coef);
        elements.emplace_back(row, idx(j-2*ed), -7.5*coef);
        elements.emplace_back(row, idx(j-3*ed), 1.5*coef);
    } else {
        std::cerr << "[Error] addElementForHomoDirichlet:: out of range!" << std::endl;
        exit(-1);
    }
}

void GePUP::initialize(){
    // Initialize the AMG solvers
    std::vector<Triple> elements;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            int row = idx(i,j);
            idpair id(i,j);
            addElementForHomoDirichlet(elements, row, id, -5.0/(dH*dH));
            for(int d = 0; d < 2; d++){
                idpair ed; ed[d]=1;
                addElementForHomoDirichlet(elements, row, id-2*ed, -1.0/(12.0*dH*dH));
                addElementForHomoDirichlet(elements, row, id-ed, 4.0/(3.0*dH*dH));
                addElementForHomoDirichlet(elements, row, id+ed, 4.0/(3.0*dH*dH));
                addElementForHomoDirichlet(elements, row, id+2*ed, -1.0/(12.0*dH*dH));
            }
        }
    poissonDirichlet.generateGrid(SparseMatrix(M*M, M*M, elements));
    
    elements.clear();
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            int row = idx(i,j);
            idpair id(i,j);
            addElementForNeumann(elements, row, id, -5.0/(dH*dH));
            for(int d = 0; d < 2; d++){
                idpair ed; ed[d]=1;
                addElementForNeumann(elements, row, id-2*ed, -1.0/(12.0*dH*dH));
                addElementForNeumann(elements, row, id-ed, 4.0/(3.0*dH*dH));
                addElementForNeumann(elements, row, id+ed, 4.0/(3.0*dH*dH));
                addElementForNeumann(elements, row, id+2*ed, -1.0/(12.0*dH*dH));
            }
        }
    poissonNeumann.generateGrid(SparseMatrix(M*M, M*M, elements));
    poissonNeumann.setPureNeumann();
}