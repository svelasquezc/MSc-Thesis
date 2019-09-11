#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include <string>

namespace Initial_Conditions{

    std::string timestamp="";
    double mytime=0;
    double simulationtime = 86400;
    double timedelta=1;
    int wells_quantity=0;
    int term=0;
    int fluids_quantity=0;
    int stencil[2] = {-1,1};
    int equilibrium_relations_quantity=0;
    int interfluid_interactions_quantity=0;
    int cells_number=0;
    int changing_wells=0;
};

#endif /* INITIAL_CONDITIONS_H */
