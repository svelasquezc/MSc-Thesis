#ifndef CELL_H
#define CELL_H

#include<vector>
#include<memory>

#include "Face.h"

class Cell{

    using Faces_t = std::vector<std::shared_ptr<Face>>;
    
 private:
    int _index;
    
    Faces_t _active_faces;
    std::vector<int> _numeration_3d;
    double _volume;
    double _depth;
    int _active_faces_quantity;
 public:

    using Face_iterator = Faces_t::iterator;
    using Face_const_iterator = Faces_t::const_iterator;

    Cell(int index): _index(index){
        _numeration_3d=std::vector<int>(3);
        _active_faces_quantity=0;
        _active_faces = std::vector<std::shared_ptr<Face>>();
    };
    
    void volume(double volume){_volume=volume;};
    void volume(double dz, double dy, double dx){_volume=dx*dy*dz;};
    const double volume() const {return _volume;};
    
    void numeration3D(int pos, int val){_numeration_3d[pos]=val;};
    const std::vector<int>& numeration3D () const {return _numeration_3d;};

    void index(int index){_index=index;};
    const int& index() const{return _index;};
    
    void depth(double depth){_depth=depth;};
    const double& depth() const{return _depth;};

    void activeFacesQuantity(int activeFacesQuantity){_active_faces_quantity=activeFacesQuantity;};
    const int& activeFacesQuantity() const{return _active_faces_quantity;};

    
    void pushFace(std::shared_ptr<Face>& face){
        _active_faces.push_back(face);
    };

    Face_iterator begin() {return _active_faces.begin();};
    Face_iterator end()   {return _active_faces.end();};

    Face_const_iterator begin()  const {return _active_faces.begin();};
    Face_const_iterator end()    const {return _active_faces.end();};
    Face_const_iterator cbegin() const {return _active_faces.cbegin();};
    Face_const_iterator cend()   const {return _active_faces.cend();};
    
};

void Face::neighborCell(std::shared_ptr<Cell>& cell) {_neighbor_cell = cell;};

#endif /* CELL_H */
