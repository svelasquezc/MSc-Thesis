#ifndef CELL_H
#define CELL_H

#include<vector>

#include "Face.h"

class Cell{
 private:
    int index;
    
    std::vector<Face> active_faces;
    std::vector<int> numeration3D;
    double Volume;
    double depth;
 public:
    int Number_of_active_faces;
    Cell(int _index): index(_index){
        numeration3D=std::vector<int>(3);
        Number_of_active_faces=0;
        active_faces = std::vector<Face>();
    };
    void setVolume(double _Volume){Volume=_Volume;};
    void setVolume(double dz, double dy, double dx){Volume=dx*dy*dz;};
    double getVolume(){return Volume;};
    
    void setNumeration3D(int pos, int val){numeration3D[pos]=val;};
    std::vector<int> getNumeration3D(){return numeration3D;};
    
    void setDepth(double _depth){depth=_depth;};
    double getDepth(){return depth;};
    void pushFace(Face _face){
        active_faces.push_back(_face);
    };
};

void Face::setNeighbor(Cell _cell){
    neighbor_cell = std::make_shared<Cell> (_cell);
};

#endif /* CELL_H */
