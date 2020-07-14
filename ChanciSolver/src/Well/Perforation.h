#ifndef PERFORATION_H
#define PERFORATION_H

#include "Constants.h"

class Perforation{

 protected:
    
    int _index;
    int _local_index;
    std::vector<int> _location = std::vector<int>(3);
    std::vector<double> _well_index;
    double _skin;
    
    std::vector<double> _equivalent_radius;

    Perforation(){
        _well_index = std::vector<double>(2,0.0);
        _equivalent_radius = std::vector<double>(2,0.0);
    };
    
 public:
    //pure virtual function: Makes perforation an abstract class
    virtual const double totalFlow() const = 0;
    virtual const std::string type() const = 0;
    
    void index(int index){_index = index;};
    void localIndex(int local_index){_local_index = local_index;};
    void location(const int x, const int y, const int z){
        _location[0]=x;
        _location[1]=y;
        _location[2]=z;
    };
    void skin(double skin){_skin = skin;};

    const int& index() const {return _index;};
    const int& local_index() const {return _local_index;};
    const double& skin() const {return _skin;};
    const int& location(const int direction) const {return _location[direction];};
    const double& equivalentRadius(const int& term) const {return _equivalent_radius[term];};
    const double& wellIndex(const int& term) const {return _well_index[term];};

    void calculateEquivalentRadius(const int& term, const double y_direction_absolute_permeability, const double x_direction_absolute_permeability, const double y_thickness, const double x_thickness){
        _equivalent_radius[term] = 0.14*sqrt(sqrt(y_direction_absolute_permeability/x_direction_absolute_permeability)*pow(x_thickness, 2.0) +
                                             std::sqrt(x_direction_absolute_permeability/y_direction_absolute_permeability) *
                                             std::pow(y_thickness, 2.0))/(0.5*(std::pow(y_direction_absolute_permeability/x_direction_absolute_permeability,1.0/4.0) +
                                                                               std::pow(x_direction_absolute_permeability/y_direction_absolute_permeability,1.0/4.0)));
    };

    void calculateWellIndex(const int& term, const double z_thickness, const double x_direction_absolute_permeability, const double y_direction_absolute_permeability, const double well_radius){
        _well_index[term]=2*pi()*z_thickness*std::sqrt(x_direction_absolute_permeability*y_direction_absolute_permeability)/(_skin+std::log(_equivalent_radius[term]/well_radius));
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



#endif /* PERFORATION_H */
