#ifndef ROCK_H
#define ROCK_H

#include "Value_Reader.h"

class Rock : protected Value_Reader{
 private:
    
    double reference_pressure;
    double compressibility;
    std::vector<std::vector<std::vector<double>>> absolute_permeability;   
    std::vector<std::vector<double>> porosity;

 public:

    Rock(){};
    void characterize(int& _cells_number);
    void calculateLinearCompressibility(int& _term, int& _cells_number, double pressure);
    void updateProperties(int& _term);
};

void Rock::characterize(int& _cells_number){
    
    std::ostringstream ss = std::ostringstream();
    const std::string axisnames[3]={"x", "y", "z"};
    
    absolute_permeability = std::vector<std::vector<std::vector<double>>>
        (1,std::vector<std::vector<double>>(_cells_number,std::vector<double>(3)));
    porosity              = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));

    myRead(std::string("Please insert rock compressibility [1/Pa]"), compressibility, std::string("Please insert a valid input"));

    myRead(std::string("Please insert reference pressure [Pa]"), reference_pressure, std::string("Please insert a valid input"));

    for(int cellindex=0; cellindex<_cells_number; ++cellindex){
        
        ss << "Please insert initial porosity for the "<< cellindex+1 << " cell [-]";
        myRead(ss.str(), porosity[0][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        
        
    };
    for(int cellindex=0; cellindex<_cells_number; ++cellindex){
        for(int direction=0; direction<3;++direction){
            ss << "Please insert initial absolute permeability for the "<< cellindex+1
               << " cell in direction " << axisnames[direction] << " [m2]";
            myRead(ss.str(), absolute_permeability[0][cellindex][direction], std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
        };
    };
}

void updateProperties(int& _term){

};

void Rock::calculateLinearCompressibility(int& _term, int& _cell_index, double pressure){
    porosity[_term][_cell_index] = porosity[0][_cell_index]
        * ( 1 + compressibility * (pressure - reference_pressure));
};

#endif /* ROCK_H */

