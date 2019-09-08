#ifndef ROCK_H
#define ROCK_H


#include <fstream>
#include <string>
#include <map>
#include <algorithm>

#include "Value_Reader.h"

class Rock : protected Value_Reader{
 private:
    
    double _reference_pressure;
    double _compressibility;
    std::map<const std::string, std::vector<std::vector<double>>> _absolute_permeability;   
    std::map<const std::string, std::vector<double>> _porosity;

 public:

    Rock(){};
    void characterize(const int& cells_number);
    void characterizeFromFile(std::ifstream& rock_reader, const int& cells_number);
    void porosity(const std::string& term, const int cell_index, const double pressure);
    
    void updateProperties(const std::string& term);

    const double& porosity (const std::string& term, const int cells_number) const {
        return _porosity.find(term)->second[cells_number];
    };

    const std::vector<double>& absolutePermeability(const std::string& term, const int cells_number) const {
        return _absolute_permeability.find(term)->second[cells_number];
    };
};

void Rock::characterize(const int& cells_number){
    
    std::ostringstream ss = std::ostringstream();
    const std::string axisnames[3]={"x", "y", "z"};
    
    _absolute_permeability["N"] = std::vector<std::vector<double>>(cells_number,std::vector<double>(3));
    _absolute_permeability["K"] = std::vector<std::vector<double>>(cells_number,std::vector<double>(3));
    _porosity["0"]              = std::vector<double>(cells_number);
    _porosity["N"]              = std::vector<double>(cells_number);
    _porosity["K"]              = std::vector<double>(cells_number);

    myRead(std::string("Please insert rock compressibility [1/Pa]"), _compressibility, std::string("Please insert a valid input"));

    myRead(std::string("Please insert reference pressure [Pa]"), _reference_pressure, std::string("Please insert a valid input"));

    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        
        ss << "Please insert initial porosity for the "<< cellindex+1 << " cell [-]";
        myRead(ss.str(), _porosity["0"][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        
        
    };
    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        for(int direction=0; direction<3;++direction){
            ss << "Please insert initial absolute permeability for the "<< cellindex+1
               << " cell in direction " << axisnames[direction] << " [m2]";
            myRead(ss.str(), _absolute_permeability["K"][cellindex][direction], std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
        };
    };
};

void Rock::characterizeFromFile(std::ifstream& rock_reader, const int& cells_number){
    std::string element;
    
    _absolute_permeability["N"] = std::vector<std::vector<double>>(cells_number,std::vector<double>(3));
    _absolute_permeability["K"] = std::vector<std::vector<double>>(cells_number,std::vector<double>(3));
    _porosity["0"]              = std::vector<double>(cells_number);
    _porosity["N"]              = std::vector<double>(cells_number);
    _porosity["K"]              = std::vector<double>(cells_number);
    
    while(rock_reader>>element){
        std::transform(element.begin(), element.end(),element.begin(), ::toupper);
        if(element == "COMPRESSIBILITY"){
            rock_reader>>_compressibility;
        }else if(element == "REFERENCE_PRESSURE"){
            rock_reader>>_reference_pressure;
        }else if(element == "POROSITY"){
            for(int cellindex=0; cellindex<cells_number; ++cellindex){
                rock_reader>>_porosity["0"][cellindex];
            };
        }else if(element == "ABSOLUTE_PERMEABILITY"){
            for(int cellindex=0; cellindex<cells_number; ++cellindex){
                for(int direction=0; direction<3;++direction){
                    rock_reader>>_absolute_permeability["N"][cellindex][direction];
                };
            };
            return;
        };
    };
};

void Rock::updateProperties(const std::string& term){
    if(term == "K"){
        _absolute_permeability["K"]=_absolute_permeability["N"];
        _porosity["K"]=_porosity["N"];
    }else{
        _absolute_permeability["N"]=_absolute_permeability["K"];
        _porosity["N"]=_porosity["K"];
    }
};

void Rock::porosity(const std::string& term, const int cell_index, const double pressure){
    _porosity[term][cell_index] = _porosity["0"][cell_index];
    //* ( 1.0 + (_compressibility * (pressure - _reference_pressure)));
};

#endif /* ROCK_H */

