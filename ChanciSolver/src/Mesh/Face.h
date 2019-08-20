#ifndef FACE_H
#define FACE_H

#include <vector>
#include <memory>

class Cell;

class Face{
 private:
    int _index;
    std::shared_ptr<Cell> _neighbor_cell;
    double _area;
    int _orientation;
 public:

    void area(double area) {_area=area;};
    const double& area() const {return _area;};
    
    void orientation(double orientation) {_orientation=orientation;};
    const int& orientation() const {return _orientation;};
    
    void index(int index){_index = index;};
    const int& index() const {return _index;};
    
    void neighbor(std::shared_ptr<Cell>& _cell);
    const std::shared_ptr<Cell>& neighbor() const {return _neighbor_cell;};
};

#endif /* FACE_H */
