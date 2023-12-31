#ifndef MESH_H
#define MESH_H

#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include <exception>
#include <string>
#include <fstream>
#include <algorithm>

#include "Cell.h"
#include "Value_Reader.h"

class Mesh : protected Value_Reader{
    
 private:
    
    using Cells_t = std::vector<std::shared_ptr<Cell>>;
    int _dimension;
    std::vector<int> _cells_quantity = std::vector<int>(3);
    int _cell_total;
    std::vector<std::vector<double>> _thickness;
    std::vector<std::vector<double>> _top;
    int _defined=0;
    Cells_t _cells;
    
 public:
    
    using Cell_iterator = Cells_t::iterator;
    using Cell_const_iterator = Cells_t::const_iterator;
    
    Mesh();
    int getCellTotal();
    void define();
    void defineFromFile(std::ifstream& mesh_reader);
    void appear(const std::string& _timestamp, const int stencil[2]);
    //int listCell(int index);
    int listCell(int posx, int posy, int posz);
    int listCell(std::vector<int> _Numeration);

    const std::shared_ptr<Cell>& cell(const int index) const {return _cells[index];};

    const double& thickness(const int axis, int spacing) const {
        if(spacing == _cells_quantity[axis]){
            --spacing;
        };
        return _thickness[axis][spacing];
    };

    const double& top(int y_index, int x_index) const {
        if(x_index == _cells_quantity[0]){
            --x_index;
        };
        if(y_index == _cells_quantity[1]){
            --y_index;
        };
        return _top[y_index][x_index];
    };
    
    const std::vector<int>& cellsQuantity() const {return _cells_quantity;};

    Cell_iterator begin() {return _cells.begin();};
    Cell_iterator end()   {return _cells.end();};

    Cell_const_iterator begin()  const {return _cells.begin();};
    Cell_const_iterator end()    const {return _cells.end();};
    Cell_const_iterator cbegin() const {return _cells.cbegin();};
    Cell_const_iterator cend()   const {return _cells.cend();};

    //const Cells_t& cells() const {return _cells;};
};

Mesh::Mesh(){};

int Mesh::getCellTotal(){return _cell_total;};

void Mesh::define(){

    const std::string axisnames[3]={"x", "y", "z"};
    std::ostringstream ss = std::ostringstream();
    int axis=1;
    double aux_thickness;
    _cell_total=1;
    _thickness = std::vector<std::vector<double>>(3,std::vector<double>());
    myRead(std::string("Please insert the dimension of the mesh: "), _dimension, std::string("Please insert a valid input"));
    for(axis=0; axis<3; ++axis){
        ss << "Please insert number of cells in direction " << axisnames[axis] << "(integer): ";
        myRead<>(ss.str(), _cells_quantity[axis], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        for(int spacing=0; spacing<_cells_quantity[axis];++spacing){
            ss << "Please insert "<< spacing << " spacing in direction " << axisnames[axis];
            myRead<>(ss.str(), aux_thickness, std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
            _thickness[axis].push_back(aux_thickness);
        };
        _cell_total=_cell_total*_cells_quantity[axis];
    };
    
    _top = std::vector<std::vector<double>>(_cells_quantity[1], std::vector<double>(_cells_quantity[0]));
    
    for(int y_spacing=0; y_spacing<_cells_quantity[1];++y_spacing){
	
	for(int x_spacing=0; x_spacing<_cells_quantity[0];++x_spacing){
	    ss << "Please insert mesh top in position (x,y): (" << x_spacing << " , " << y_spacing << ") (double): ";
	    myRead<>(ss.str(), _top[y_spacing][x_spacing], std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
	};
	
    };
    _defined = 1;
};

void Mesh::defineFromFile(std::ifstream& mesh_reader){

    std::string element;
    int axis=1;
    double aux_thickness;

    _cell_total=1;
    _thickness = std::vector<std::vector<double>>(3,std::vector<double>());

    while(mesh_reader>> element){
        std::transform(element.begin(), element.end(),element.begin(), ::toupper);
        if(element == "DIMENSIONS"){
            for(axis=0; axis<3; ++axis){
                mesh_reader >>_cells_quantity[axis];
                _cell_total=_cell_total*_cells_quantity[axis];
            };
        }else if(element == "THICKNESS"){
            for(axis=0; axis<3; ++axis){
                _thickness[axis].reserve(_cells_quantity[axis]);
                for(int spacing=0; spacing<_cells_quantity[axis];++spacing){
                    mesh_reader >> aux_thickness;
                    _thickness[axis].push_back(aux_thickness);
                };
            };
        }else if(element == "TOPS"){
            _top = std::vector<std::vector<double>>(_cells_quantity[1], std::vector<double>(_cells_quantity[0]));
            for(int y_spacing=0; y_spacing<_cells_quantity[1];++y_spacing){
                for(int x_spacing=0; x_spacing<_cells_quantity[0];++x_spacing){
                    mesh_reader >> _top[y_spacing][x_spacing];
                };
            };
            break;
        };
    };
    _defined = 1;
    return;
};

void Mesh::appear(const std::string& _timestamp, const int stencil[2]){
    if(_timestamp=="" && _defined==1){
        std::shared_ptr<Cell> my_cell;
        int index = 0;
        double auxiliary_centroid=0;
        for(int z_position=0;z_position<_cells_quantity[2];++z_position){
	    
            if(z_position==0){
                auxiliary_centroid += _thickness[2][z_position]/2.0;
            }else{
                auxiliary_centroid += _thickness[2][z_position-1]/2.0 + _thickness[2][z_position]/2.0;
            };
	    
            for(int y_position=0;y_position<_cells_quantity[1];++y_position){
                for(int x_position=0;x_position<_cells_quantity[0];++x_position){
                    my_cell = std::make_shared<Cell>(index);
                    my_cell->volume(_thickness[2][z_position], _thickness[1][y_position], _thickness[0][x_position]);
                    my_cell->depth(auxiliary_centroid + _top[y_position][x_position]);
                    my_cell->numeration3D(0,x_position);
                    my_cell->numeration3D(1,y_position);
                    my_cell->numeration3D(2,z_position);
                    _cells.push_back(std::move(my_cell));
                    ++index;
                };
            };
        };
        int face_index=0;
        for(auto celli : _cells){
            
            auto local_numeration = celli->numeration3D();
            
            for(int axis=0; axis<3; ++axis){
                for(int direction=0;direction<2;++direction){
                    local_numeration[axis]=local_numeration[axis]+stencil[direction];
                    int neighbor_cell=listCell(local_numeration);
                    
                    local_numeration[axis]=local_numeration[axis]-stencil[direction];
                    
                    if(neighbor_cell > -1){
                        
                        celli->activeFacesQuantity(celli->activeFacesQuantity()+1);
                        std::shared_ptr<Face> face = std::make_shared<Face>();
                        switch(axis){
                        case 0:
                            face->area(_thickness[1][local_numeration[1]]*_thickness[2][local_numeration[2]]);
                            face->orientation(0);
                            break;
                        case 1:
                            face->area(_thickness[0][local_numeration[0]]*_thickness[2][local_numeration[2]]);
                            face->orientation(1);
                            break;
                        case 2:
                            face->area(_thickness[0][local_numeration[0]]*_thickness[1][local_numeration[1]]);
                            face->orientation(2);
                            break;
                        }
                        face->neighborCell(_cells[neighbor_cell]);
                        ++face_index;
                        face->index(face_index);
                        celli->pushFace(face);
                    }
                    
                }
            }
        }
        _defined = 2;
    };
};

int Mesh::listCell(int posx, int posy, int posz){
    int index;
    if (!(posx < 0 || posy < 0 || posz < 0 || posx >= _cells_quantity[0] || posy >= _cells_quantity[1] || posz >= _cells_quantity[2])){
        index = posx + posy*_cells_quantity[0] + posz*_cells_quantity[0]*_cells_quantity[1];
        return index;
    }else{
        return -1;
    }
}

int Mesh::listCell(std::vector<int> numeration){
    int index;
    if (!(numeration[0] < 0 || numeration[1] < 0 || numeration[2] < 0 || numeration[0] >= _cells_quantity[0] || numeration[1] >= _cells_quantity[1] || numeration[2] >= _cells_quantity[2])){
        index = numeration[0] + numeration[1]*_cells_quantity[0] + numeration[2]*_cells_quantity[0]*_cells_quantity[1];
        return index;
    }else{
        return -1;
    }
    
}

#endif /* MESH_H */
