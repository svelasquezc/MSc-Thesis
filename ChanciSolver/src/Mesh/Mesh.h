#ifndef MESH_H
#define MESH_H

#include <vector>
#include <iostream>
#include <sstream>
#include <exception>
#include <string>

#include "Cell.h"
#include "Value_Reader.h"

class Mesh : protected Value_Reader{
    
 private:
    
    using Cells_t = std::vector<Cell>;
    int _dimension;
    std::vector<int> _cell_number = std::vector<int>(3);
    int _cell_total;
    std::vector<std::vector<double>> _thickness;
    mutable bool _defined=false;
    Cells_t _cells;
    
 public:
    
    using Cell_iterator = Cells_t::iterator;
    using Cell_const_iterator = Cells_t::const_iterator;
    
    Mesh();
    int getCellTotal();
    void defineMesh();
    void appear(const std::string& _timestamp, const int stencil[2]);
    //int listCell(int index);
    int listCell(int posx, int posy, int posz);
    int listCell(std::vector<int> _Numeration);

    const double& thickness(const int axis, const int spacing) const {return _thickness[axis][spacing];};

    Cell_iterator begin() {return _cells.begin();};
    Cell_iterator end()   {return _cells.end();};

    Cell_const_iterator begin()  const {return _cells.begin();};
    Cell_const_iterator end()    const {return _cells.end();};
    Cell_const_iterator cbegin() const {return _cells.cbegin();};
    Cell_const_iterator cend()   const {return _cells.cend();};
};

Mesh::Mesh(){};

int Mesh::getCellTotal(){return _cell_total;};

void Mesh::defineMesh(){

    const std::string axisnames[3]={"x", "y", "z"};
    std::ostringstream ss = std::ostringstream();
    int axis=1;
    double aux_thickness;
    _cell_total=1;
    _thickness = std::vector<std::vector<double>>(3,std::vector<double>());
    myRead(std::string("Please insert the dimension of the mesh: "), _dimension, std::string("Please insert a valid input"));
    for(axis=0; axis<3; axis++){
        ss << "Please insert number of cells in direction " << axisnames[axis] << "(integer): ";
        myRead<>(ss.str(), _cell_number[axis], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        for(int spacing=0; spacing<_cell_number[axis];spacing++){
            ss << "Please insert "<< spacing << " spacing in direction " << axisnames[axis];
            myRead<>(ss.str(), aux_thickness, std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
            _thickness[axis].push_back(aux_thickness);
        };
        _cell_total=_cell_total*_cell_number[axis];
    };
    _defined = true;
};

void Mesh::appear(const std::string& _timestamp, const int stencil[2]){
    if(_timestamp=="" && _defined){
        int index = 0;
        for(int axisz=0;axisz<_cell_number[2];axisz++){
            for(int axisy=0;axisy<_cell_number[1];axisy++){
                for(int axisx=0;axisx<_cell_number[0];axisx++){
                    Cell my_cell = Cell(index);
                    my_cell.volume(_thickness[2][axisz], _thickness[1][axisy], _thickness[0][axisx]);
                    my_cell.numeration3D(0,axisx);
                    my_cell.numeration3D(1,axisy);
                    my_cell.numeration3D(2,axisz);
                    _cells.push_back(my_cell);
                    index++;
                };
            };
        };
        int face_index=0;
        for(auto celli : _cells){
            
            auto local_numeration = celli.numeration3D();
            
            for(int axis=0; axis<3; axis++){
                for(int direction=0;direction<2;direction++){
                    local_numeration[axis]=local_numeration[axis]+stencil[direction];
                    int neighbor_cell=listCell(local_numeration);
                    
                    local_numeration[axis]=local_numeration[axis]-stencil[direction];
                    
                    if(neighbor_cell > -1){
                        
                        celli.numberOfActiveFaces(celli.numberOfActiveFaces()+1);
                        Face face = Face();
                        switch(axis){
                        case 0:
                            face.area(_thickness[1][local_numeration[1]]*_thickness[2][local_numeration[2]]);
                            face.orientation(0);
                            break;
                        case 1:
                            face.area(_thickness[0][local_numeration[0]]*_thickness[2][local_numeration[2]]);
                            face.orientation(1);
                            break;
                        case 2:
                            face.area(_thickness[0][local_numeration[0]]*_thickness[1][local_numeration[1]]);
                            face.orientation(2);
                            break;
                        }
                        face.neighbor(_cells[neighbor_cell]);
                        face_index++;
                        face.index(face_index);
                        celli.pushFace(face);
                    }
                    
                }
            }
        }
        
    };
};

int Mesh::listCell(int posx, int posy, int posz){
    int index;
    if (!(posx < 0 || posy < 0 && posz < 0 || posx >= _cell_number[0] || posy >= _cell_number[1] || posz >= _cell_number[2])){
        index = posx + posy*_cell_number[0] + posz*_cell_number[0]*_cell_number[1];
	return index;
    }else{
        return -1;
    }
}

int Mesh::listCell(std::vector<int> numeration){
    int index;
    if (!(numeration[0] < 0 || numeration[1] < 0 || numeration[2] > 0 || numeration[0] >= _cell_number[0] || numeration[1] >= _cell_number[1] || numeration[2] >= _cell_number[2])){
        index = numeration[0] + numeration[1]*_cell_number[0] + numeration[2]*_cell_number[0]*_cell_number[1];
	return index;
    }else{
        return -1;
    }
    
}

#endif /* MESH_H */
