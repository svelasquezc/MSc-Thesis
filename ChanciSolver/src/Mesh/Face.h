#ifndef FACE_H
#define FACE_H

#include <vector>
#include <memory>

class Cell;

class Face{
 private:
    int index;
    std::shared_ptr<Cell> neighbor_cell;
    double area;
    int orientation;
 public:

    void setArea(double _area) {area=_area;};
    void setOrientation(double _orientation) {orientation=_orientation;};
    void setIndex(int _index){index = _index;};
    double getArea(){return area;};
    int getOrientation(){return orientation;};
    void setNeighbor(Cell _cell);    
};

#endif /* FACE_H */
