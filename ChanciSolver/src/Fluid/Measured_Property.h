#ifndef MEASURED_PROPERTY_H
#define MEASURED_PROPERTY_H

#include <map>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "Value_Reader.h"

class Measured_Property : protected Value_Reader{
 private:

    using Measures_t = std::map<double,double>;
    
    
    std::string _reference_measure_name;
    std::string _measured_value_name;
    int _number_of_measures;
    Measures_t _taken_measures;
    
 public:
    
    using Measures_iterator = Measures_t::iterator;
    using Measures_const_iterator = Measures_t::const_iterator;
    
 Measured_Property(std::string reference_measure_name, std::string measured_value_name):
    _reference_measure_name(reference_measure_name), _measured_value_name(measured_value_name)
    {
        _taken_measures = std::map<double,double>();
    };
    void readMe();
    double interpolate(double reference_value);

    Measures_iterator begin() {return _taken_measures.begin();};
    Measures_iterator end()   {return _taken_measures.end();};

    Measures_const_iterator begin()  const {return _taken_measures.begin();};
    Measures_const_iterator end()    const {return _taken_measures.end();};
    Measures_const_iterator cbegin() const {return _taken_measures.cbegin();};
    Measures_const_iterator cend()   const {return _taken_measures.cend();};
};

void Measured_Property::readMe(){
    std::ostringstream ss=std::ostringstream();
    double aux_reference, aux_measure;
    
    ss << "Please insert the number of measures for " << _measured_value_name <<": ";
    myRead(ss.str(), _number_of_measures,std::string("Please insert a valid input"));
    ss.str("");
    ss.clear();
    
    for(int measure=1; measure<=_number_of_measures; ++measure){
	ss<<"Please insert the " << measure << " " << _reference_measure_name;
	myRead(ss.str(), aux_reference,std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();

	ss<<"Please insert the " << measure << " " << _measured_value_name;
	myRead(ss.str(), aux_measure,std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();

	_taken_measures.insert(std::make_pair(aux_reference,aux_measure));
    };

};
/*
  This Function implements Linear interpolation in the taken_measures map
*/
double Measured_Property::interpolate(double reference_value){

    auto great_it = _taken_measures.upper_bound(reference_value);

    if(great_it == _taken_measures.end()){
        return (--great_it)->second;
    }
    if(great_it == _taken_measures.begin()){
        return great_it->second;
    }

    auto less_it = great_it;
    --less_it;
    
    //Interpolate
    return less_it->second + (reference_value - less_it->first) * (great_it->second - less_it->second) / (great_it->first - less_it->first);
};

#endif /* MEASURED_PROPERTY_H */
