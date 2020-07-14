#ifndef PROPERTY_MEASURE_H
#define PROPERTY_MEASURE_H

#include <map>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>

#include "Value_Reader.h"

class Property_Measure : protected Value_Reader{
private:

    using Measures_t = std::map<double,double>;
    
    
    std::string _reference_measure_name;
    std::string _measured_value_name;
    int _measures_quantity;
    Measures_t _taken_measures;
    
public:
    
    using Measures_iterator = Measures_t::iterator;
    using Measures_const_iterator = Measures_t::const_iterator;
    
    Property_Measure(std::string reference_measure_name, std::string measured_value_name):
        _reference_measure_name(reference_measure_name), _measured_value_name(measured_value_name)
    {
        _taken_measures = std::map<double,double>();
    };
    void readMe();
    void readFromFile(std::ifstream& property_reader);
    double interpolate(double reference_value);

    Measures_iterator begin() {return _taken_measures.begin();};
    Measures_iterator end()   {return _taken_measures.end();};

    Measures_const_iterator begin()  const {return _taken_measures.begin();};
    Measures_const_iterator end()    const {return _taken_measures.end();};
    Measures_const_iterator cbegin() const {return _taken_measures.cbegin();};
    Measures_const_iterator cend()   const {return _taken_measures.cend();};
};

void Property_Measure::readMe(){
    std::ostringstream ss=std::ostringstream();
    double aux_reference, aux_measure;
    
    ss << "Please insert the number of measures for " << _measured_value_name <<": ";
    myRead(ss.str(), _measures_quantity,std::string("Please insert a valid input"));
    ss.str("");
    ss.clear();
    
    for(int measure=1; measure<=_measures_quantity; ++measure){
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

void Property_Measure::readFromFile(std::ifstream& property_reader){

    double reference_measure, value_measure;
    
    property_reader>> _measures_quantity;

    for (int measure=0; measure<_measures_quantity; ++measure){
        property_reader>>reference_measure;
        property_reader>>value_measure;

        _taken_measures.insert(std::make_pair(reference_measure,value_measure));
    };
    return;
};

/*
  This Function implements Linear interpolation in the taken_measures map
*/
double Property_Measure::interpolate(double reference_value){

    auto great_it = _taken_measures.lower_bound(reference_value);

    if(great_it == _taken_measures.end()){
        --great_it;
        return great_it->second;
    }
    if(great_it == _taken_measures.begin()){
        return great_it->second;
    }

    auto less_it = great_it;
    --less_it;
    
    //Interpolate

    auto m = (great_it->second - less_it->second) / (great_it->first - less_it->first);

    auto b = less_it->second - m*less_it->first;
    
    return m*reference_value + b;
};

#endif /* PROPERTY_MEASURE_H */
