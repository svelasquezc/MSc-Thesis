#ifndef ROCK_H
#define ROCK_H

#include "Value_Reader.h"

class Rock : protected Value_Reader{
 private:
    
    double _reference_pressure;
    double _compressibility;
    std::vector<std::vector<std::vector<double>>> _absolute_permeability;   
    std::vector<std::vector<double>> _porosity;

 public:

    Rock(){};
    void characterize(int& _cells_number);
    void porosity(const int term, const int cell_index, const double pressure);
    void updateProperties(int& _term);

    const double& porosity (const int term, const int cells_number) const {return _porosity[term][cells_number];};
};

void Rock::characterize(int& cells_number){
    
    std::ostringstream ss = std::ostringstream();
    const std::string axisnames[3]={"x", "y", "z"};
    
    _absolute_permeability = std::vector<std::vector<std::vector<double>>>
        (1,std::vector<std::vector<double>>(cells_number,std::vector<double>(3)));
    _porosity              = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));

    myRead(std::string("Please insert rock compressibility [1/Pa]"), _compressibility, std::string("Please insert a valid input"));

    myRead(std::string("Please insert reference pressure [Pa]"), _reference_pressure, std::string("Please insert a valid input"));

    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        
        ss << "Please insert initial porosity for the "<< cellindex+1 << " cell [-]";
        myRead(ss.str(), _porosity[0][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        
        
    };
    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        for(int direction=0; direction<3;++direction){
            ss << "Please insert initial absolute permeability for the "<< cellindex+1
               << " cell in direction " << axisnames[direction] << " [m2]";
            myRead(ss.str(), _absolute_permeability[0][cellindex][direction], std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
        };
    };
}

void updateProperties(int& term){

};

void Rock::porosity(const int term, const int cell_index, const double pressure){
    _porosity[term][cell_index] = _porosity[0][cell_index]
        * ( 1 + _compressibility * (pressure - _reference_pressure));
};

#endif /* ROCK_H */

