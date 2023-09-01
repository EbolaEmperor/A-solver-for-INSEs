#include "GePUP_IMEX.h"

namespace ERK_ESDIRK_Table{
    const int stage = 6;
    const double gma = 0.25;
    const double b[6] = {
        0.15791629516167136, 0, 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, gma
    };
    const double aE[6][6] = {
        0, 0, 0, 0, 0, 0,
        2.0*gma, 0, 0, 0, 0, 0,
        0.221776, 0.110224, 0, 0, 0, 0,
        -0.04884659515311857, -0.17772065232640102, 0.8465672474795197, 0, 0, 0,
        -0.15541685842491548, -0.3567050098221991, 1.0587258798684427, 0.30339598837867193, 0, 0,
        0.2014243506726763, 0.008742057842904185, 0.15993995707168115, 0.4038290605220775, 0.22606457389066084, 0
    };
    const double aI[6][6] = {
        0, 0, 0, 0, 0, 0,
        gma, gma, 0, 0, 0, 0,
        0.137776, -0.055776, gma, 0, 0, 0,
        0.14463686602698217, -0.22393190761334475, 0.4492950415863626, gma, 0, 0,
        0.09825878328356477, -0.5915442428196704, 0.8101210538282996, 0.283164405707806, gma, 0,
        b[0], b[1], b[2], b[3], b[4], b[5]
    };
    const double c[6] = {
        0, 0.5, 0.332, 0.62, 0.85, 1.0
    };
}

void GePUP_IMEX::initialize(){
    GePUP::initialize();
    using namespace ERK_ESDIRK_Table;
    std::vector<Triple> elements;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            int row = idx(i,j);
            idpair id(i,j);
            addElementForHomoDirichlet(elements, row, id, 1.0+5.0*dT*nu*gma/(dH*dH));
            for(int d = 0; d < 2; d++){
                idpair ed; ed[d]=1;
                addElementForHomoDirichlet(elements, row, id-2*ed, dT*nu*gma/(12.0*dH*dH));
                addElementForHomoDirichlet(elements, row, id-ed, -dT*nu*gma*4.0/(3.0*dH*dH));
                addElementForHomoDirichlet(elements, row, id+ed, -dT*nu*gma*4.0/(3.0*dH*dH));
                addElementForHomoDirichlet(elements, row, id+2*ed, dT*nu*gma/(12.0*dH*dH));
            }
        }
    IRK.setStrongThereshold(std::min(0.2, 5*nu));
    IRK.generateGrid(SparseMatrix(M*M, M*M, elements));
}

void GePUP_IMEX::solve(){
    initialize();
    using namespace ERK_ESDIRK_Table;
    Field us[stage], ws[stage], XEus[stage];
    for(double t = 0.0; t+1e-12 < tEnd; t += dT){
        ws[0] = w;
        std::cerr << "Time: " << t << std::endl;
        XEus[0] = XE(u,t);
        for(int s = 1; s < stage; s++){
            double ts = t + c[s]*dT;
            Field rhs = w;
            for(int j = 0; j < s; j++){
                rhs += (dT*aE[s][j]) * XEus[j];
                rhs += (dT*nu*aI[s][j]) * L(ws[j]);
            }
            ws[s][0] = IRK.solve(rhs[0], "FMG", 8, eps);
            ws[s][1] = IRK.solve(rhs[1], "FMG", 8, eps);
            us[s] = Proj(ws[s]);
            XEus[s] = XE(us[s],ts);
        }
        Field wstar = ws[stage-1];
        for(int j = 0; j < stage; j++)
            wstar +=( dT*(b[j]-aE[stage-1][j]) ) * XEus[j];
        w = u = Proj(wstar);
    }
}