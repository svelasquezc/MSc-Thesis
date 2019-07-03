#ifndef EQUILIBRIUM_RELATION_H
#define EQUILIBRIUM_RELATION_H

#include "Fluid.h"

class Equilibrium_Relation{
 private:
    int index;
    std::shared_ptr<Fluid> contributor_fluid;
    std::shared_ptr<Fluid> receiver_fluid;
    std::vector<std::vector<double>> partition_coefficient;
    std::unique_ptr<Measured_Property> measured_partition_coefficient;
 public:
    Equilibrium_Relation(){};
    void add(int fluids_quantity);
};

#endif /* EQUILIBRIUM_RELATION_H */
