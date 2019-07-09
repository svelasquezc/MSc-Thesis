#ifndef MEASURED_PROPERTY_H
#define MEASURED_PROPERTY_H

#include <map>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "Value_Reader.h"

class Measured_Property : protected Value_Reader{
 private:
    std::string reference_measure_name;
    std::string measured_value_name;
    int number_of_measures;
    std::map<double,double> taken_measures;
    
 public:
 Measured_Property(std::string _reference_measure_name, std::string _measured_value_name):
    reference_measure_name(_reference_measure_name), measured_value_name(_measured_value_name)
    {
        taken_measures = std::map<double,double>();
    };
    void readMe();
    double interpolate(double Reference_value);
};

void Measured_Property::readMe(){
    std::ostringstream ss=std::ostringstream();
    double aux_reference, aux_measure;
    myRead(std::string("Please insert the number of measures: "), number_of_measures,std::string("Please insert a valid input"));
    
    for(int measure=1; measure<=number_of_measures; ++measure){
	ss<<"Please insert the " << measure << " " << reference_measure_name;
	myRead(ss.str(), aux_reference,std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();

	ss<<"Please insert the " << measure << " " << measured_value_name;
	myRead(ss.str(), aux_measure,std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();

	taken_measures.insert(std::make_pair(aux_reference,aux_measure));
    };

};
/*
  This Function implements Linear interpolation in the taken_measures map
 */
double Measured_Property::interpolate(double Reference_value){

    auto great_it = taken_measures.upper_bound(Reference_value);

    if(great_it == taken_measures.end()){
        return (--great_it)->second;
    }
    if(great_it == taken_measures.begin()){
        return great_it->second;
    }

    auto less_it = great_it;
    --less_it;
    
    //Interpolate
    return less_it->second + (Reference_value - less_it->first) * (great_it->second - less_it->second) / (great_it->first - less_it->first);
};

#endif /* MEASURED_PROPERTY_H */
