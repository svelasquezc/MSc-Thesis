#ifndef MEASURED_PROPERTY_H
#define MEASURED_PROPERTY_H

#include <map>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "Value_Reader.h"

class Measured_Property : protected Value_Reader{
 private:
    int number_of_measures;
    std::map<double,double> taken_measures;
    
 public:
    Measured_Property(){
        taken_measures = std::map<double,double>();
    };
    void readMe();
    double interpolate(double Reference_value);
};

void Measured_Property::readMe(){
    std::stringstream ss=std::stringstream();
    double auxReference, auxMeasure;
    myRead(std::string("Please insert the number of measures"), number_of_measures,std::string("Please insert a valid input"));
    
    for(int measure=1; measure<=number_of_measures; ++measure){
	ss<<"Please insert the "<< measure << " Reference Measure"<<std::endl;
	myRead(ss.str(), auxReference,std::string("Please insert a valid input"));
	ss.flush();

	ss<<"Please insert the "<< measure << " Measured Value"<<std::endl;
	myRead(ss.str(), auxMeasure,std::string("Please insert a valid input"));
	ss.flush();

	taken_measures.insert(std::make_pair(auxReference,auxMeasure));
    };

};
/*
  This Function implements Linear interpolation in the taken_measures map
 */
double Measured_Property::interpolate(double Reference_value){

    auto GreatIt = taken_measures.upper_bound(Reference_value);

    if(GreatIt == taken_measures.end()){
        return (--GreatIt)->second;
    }
    if(GreatIt == taken_measures.begin()){
        return GreatIt->second;
    }

    auto LessIt = GreatIt;
    --LessIt;
    
    //Interpolate
    return LessIt->second + (Reference_value - LessIt->first) * (GreatIt->second - LessIt->second) / (GreatIt->first - LessIt->first);
};

#endif /* MEASURED_PROPERTY_H */
