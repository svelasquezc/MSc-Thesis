#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>

namespace Global{

    std::string timestamp="";
    double mytime=0;
    double simulationtime = 60;
    double timedelta=1;
    int wells_quantity=0;
    int term=0;
    int fluids_quantity=0;
    int stencil[2] = {-1,1};
    int equilibrium_relations_quantity=0;
    int cells_number=0;
    int total_of_perforations;
};

#endif /* GLOBAL_H */
