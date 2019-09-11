#ifndef DATABASE_H
#define DATABASE_H

#include "Rock.h"
#include "Mesh.h"
#include "Well.h"
#include "Equilibrium_Relation.h"
#include "Interfluid_Interaction.h"

using namespace std;

namespace Database{
    vector<shared_ptr<Equation_Base>> equations =
        vector<shared_ptr<Equation_Base>>();

    vector<shared_ptr<Fluid>> characterized_fluids =
        vector<shared_ptr<Fluid>>();

    vector<unique_ptr<Equilibrium_Relation>> added_equilibrium_relations =
        vector<unique_ptr<Equilibrium_Relation>>();

    vector<unique_ptr<Interfluid_Interaction>> added_interfluid_interactions =
        vector<unique_ptr<Interfluid_Interaction>>();

    vector<shared_ptr<Well>> perforated_wells =
        vector<shared_ptr<Well>>();

    unique_ptr<Mesh> mymesh;
    unique_ptr<Rock> myrock;
};

#endif /* DATABASE_H */
