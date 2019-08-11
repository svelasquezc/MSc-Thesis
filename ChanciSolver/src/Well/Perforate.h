#ifndef PERFORATE_H
#define PERFORATE_H

#include <cmath>

class Perforate{

 private:
    
    int _index;
    int _local_index;
    std::vector<int> _position = std::vector<int>(3);
    std::vector<double> _well_index;
    double _skin;
    std::vector<double> _equivalent_radius;

 public:

    virtual void flow() = 0;
    
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

    void calculateEquivalentRadius(const int& term, const double permeability_in_y_direction, const double permeability_in_x_direction, const double thickness_in_y_direction, const double thickness_in_x_direction){
        double equivalent_radius = 0.14*sqrt(sqrt(permeability_in_y_direction/permeability_in_x_direction)*pow(thickness_in_x_direction, 2.0) + sqrt(permeability_in_x_direction/permeability_in_y_direction)*pow(thickness_in_y_direction, 2.0))/(0.5*(pow(permeability_in_y_direction/permeability_in_x_direction,1.0/4.0) + pow(permeability_in_x_direction/permeability_in_y_direction,1.0/4.0)));
    };

    void calculateWellIndex(){};
};



#endif /* PERFORATE_H */
