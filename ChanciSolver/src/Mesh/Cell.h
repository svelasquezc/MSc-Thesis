#ifndef CELL_H
#define CELL_H

#include<vector>

#include "Face.h"

class Cell{
 private:
    int _index;
    
    std::vector<Face> _active_faces;
    std::vector<int> _numeration3D;
    double _volume;
    double _depth;
 public:
    int _number_of_active_faces;
    Cell(int index): _index(index){
        _numeration_3d=std::vector<int>(3);
        _number_of_active_faces=0;
        _active_faces = std::vector<Face>();
    };
    
    void volume(double volume){_volume=volume;};
    void volume(double dz, double dy, double dx){_volume=dx*dy*dz;};
    const double volume() const {return _volume;};
    
    void numeration3D(int pos, int val){_numeration3D[pos]=val;};
    const std::vector<int>& numeration3D const(){return _numeration3D;};
    
    void depth(double depth){_depth=depth;};
    const double& depth() const{return _depth;};
    
    void pushFace(Face face){
        _active_faces.push_back(face);
    };
};

void Face::setNeighbor(Cell cell){
    neighbor_cell = std::make_shared<Cell> (cell);
};

#endif /* CELL_H */
