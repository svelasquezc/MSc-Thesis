#ifndef PERFORATE_H
#define PERFORATE_H

#include "Constants.h"

class Perforate{

 protected:
    
    int _index;
    int _local_index;
    std::vector<int> _position = std::vector<int>(3);
    std::map<const int, double> _well_index;
    double _skin;
    
    std::map<const int, double> _equivalent_radius;

    Perforate(){
        _well_index[0] = 0.0;
        _well_index[1] = 0.0;
        _equivalent_radius[0]=0.0;
        _equivalent_radius[1]=0.0;
    };
    
 public:
    //pure virtual function: Makes perforate an abstract class
    virtual const double totalFlow() const = 0;
    virtual const std::string type() const = 0;
    
    void index(int index){_index = index;};
    void localIndex(int local_index){_local_index = local_index;};
    void position(const int x, const int y, const int z){
        _position[0]=x;
        _position[1]=y;
        _position[2]=z;
    };
    void skin(double skin){_skin = skin;};

    const int& index() const {return _index;};
    const int& local_index() const {return _local_index;};
    const double& skin() const {return _skin;};
    const int& position(const int direction) const {return _position[direction];};
    const double& equivalentRadius(const int& term) const {return _equivalent_radius.find(term)->second;};
    const double& wellIndex(const int& term) const {return _well_index.find(term)->second;};

    void calculateEquivalentRadius(const int& term, const double permeability_in_y_direction, const double permeability_in_x_direction, const double thickness_in_y_direction, const double thickness_in_x_direction){
        _equivalent_radius[term] = 0.14*sqrt(sqrt(permeability_in_y_direction/permeability_in_x_direction)*pow(thickness_in_x_direction, 2.0) + std::sqrt(permeability_in_x_direction/permeability_in_y_direction)*std::pow(thickness_in_y_direction, 2.0))/(0.5*(std::pow(permeability_in_y_direction/permeability_in_x_direction,1.0/4.0) + std::pow(permeability_in_x_direction/permeability_in_y_direction,1.0/4.0)));
    };

    void calculateWellIndex(const int& term, const double z_thickness, const double x_direction_permeability, const double y_direction_permeability, const double well_radius){
        _well_index[term]=2*pi()*z_thickness*std::sqrt(x_direction_permeability*y_direction_permeability)/(_skin+std::log(_equivalent_radius[term]/well_radius));
    };

    void updateProperties(const int& term){
        if(term==1){
            _equivalent_radius[1]=_equivalent_radius[0];
            _well_index[1]=_well_index[0];
        }else{
            _equivalent_radius[0]=_equivalent_radius[1];
            _well_index[0]=_well_index[1];
        }
    };
    
};



#endif /* PERFORATE_H */
