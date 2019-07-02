#ifndef MEASURED_PROPERTY_H
#define MEASURED_PROPERTY_H

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

class Measured_Property : protected Value_Reader{
 private:
    int number_of_measures;
    std::vector<std::pair<double,double>> taken_measures;
    
 public:
    Measured_Property(){
        taken_measures = std::vector<std::pair<double,double>>();
    };
    void ReadMe();
    double interpolate(double Reference_value);
};

void Measured_Property::ReadMe(){
    std::stringstream ss=std::stringstream();
    double auxReference, auxMeasure;
    myRead(std::string("Please insert the number of measures"), &number_of_measures,std::string("Please insert a valid input"));
    
    for(auto measure : number_of_measures){
	ss<<"Please insert the "<< measure << " Reference Measure"<<std::endl;
	myRead(ss.str(), &auxReference,std::string("Please insert a valid input"));
	ss.flush();

	ss<<"Please insert the "<< measure << " Measured Value"<<std::endl;
	myRead(ss.str(), &auxMeasure,std::string("Please insert a valid input"));
	ss.flush();

	taken_measures.push_back(std::make_pair(auxReference,auxMeasure));
    };

    std::sort(taken_measures.begin(), taken_measures.end());
};
/*
  This Function implements Linear interpolation in the taken_measures vector of pairs
 */
double Measured_Property::interpolate(double Reference_value){


    auto compare = [](const std::pair<double,double>& lhs, const std::pair<double,double>& rhs){
                                      return lhs.first() < rhs.first();
    };
    
    //Find lesser adjacent
    auto lower = std::lower_bound(taken_measures.begin(), taken_measures.end(),
                                  std::pair<double,double>(Reference_value,0),
                                  compare);
    //Find greater adjacent
    auto upper = std::upper_bound(taken_measures.begin(), taken_measures.end(),
                                  std::pair<double,double>(Reference_value,0),
                                  compare);
    //Interpolate
    return lower.second() + (Reference_value-lower.first())*(upper.second()-lower.second())/(upper.first()-lower.first())
};
