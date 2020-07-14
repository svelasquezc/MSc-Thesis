#ifndef DATABASE_H
#define DATABASE_H

#include "Rock.h"
#include "Mesh.h"
#include "Well.h"
#include "Equilibrium_Relationship.h"
#include "Interphase_Interaction.h"

using namespace std;

namespace Database{
    vector<shared_ptr<Equation_Base>> equations =
        vector<shared_ptr<Equation_Base>>();

    vector<shared_ptr<Phase>> characterized_phases =
        vector<shared_ptr<Phase>>();

    vector<unique_ptr<Equilibrium_Relationship>> added_equilibrium_relationships =
        vector<unique_ptr<Equilibrium_Relationship>>();

    vector<unique_ptr<Interphase_Interaction>> added_interphase_interactions =
        vector<unique_ptr<Interphase_Interaction>>();

    vector<shared_ptr<Well>> perforated_wells =
        vector<shared_ptr<Well>>();

    unique_ptr<Mesh> defined_mesh;
    unique_ptr<Rock> characterized_rock;

    vector<double> print_times = vector<double>();
};

#endif /* DATABASE_H */
