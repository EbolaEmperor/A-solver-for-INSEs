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

ColVector GePUP::D(const TimeFunction2D *g, const double &t) const{
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) = (- g[1].intFixY_order6(j*dH, i*dH, (i+1)*dH, t)
                             + g[1].intFixY_order6((j+1)*dH, i*dH, (i+1)*dH, t)
                             - g[0].intFixX_order6(i*dH, j*dH, (j+1)*dH, t)
                             + g[0].intFixX_order6((i+1)*dH, j*dH, (j+1)*dH, t)) / (dH*dH);
        }
    return res;
}