#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>

namespace Global{

    std::string timestamp="";
    double mytime=0;
    double simulationtime = 86400;
    double timedelta=1;
    int wells_quantity=0;
    int term=1;
    int fluids_quantity=0;
    int stencil[2] = {-1,1};
    int equilibrium_relations_quantity=0;
    int interfluid_interactions_quantity=0;
    int cells_number=0;
    int changing_wells=0;
};

#endif /* GLOBAL_H */
