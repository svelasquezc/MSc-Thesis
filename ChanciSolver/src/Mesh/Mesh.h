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
    int dimension;
    std::vector<int> cell_number = std::vector<int>(3);
    int cell_total;
    std::vector<std::vector<double>> thickness;
    mutable bool defined=false;
    std::vector<Cell> cells;
    
    
 public:
    Mesh();
    int getCellTotal();
    void defineMesh();
    void appear(std::string& _timestamp, int stencil[2]);
    //int listCell(int index);
    int listCell(int posx, int posy, int posz);
    int listCell(std::vector<int> _Numeration);
};

Mesh::Mesh(){};

int Mesh::getCellTotal(){return cell_total;};

void Mesh::defineMesh(){

    const std::string axisnames[3]={"x", "y", "z"};
    std::ostringstream ss = std::ostringstream();
    int axis=1;
    double AuxThickness;
    cell_total=1;
    thickness = std::vector<std::vector<double>>(3,std::vector<double>());
    myRead(std::string("Please insert the dimension of the mesh: "), dimension, std::string("Please insert a valid input"));
    for(axis=0; axis<3; axis++){
        ss << "Please insert number of cells in direction " << axisnames[axis] << "(integer): ";
        myRead<>(ss.str(), cell_number[axis], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        for(int spacing=0; spacing<cell_number[axis];spacing++){
            ss << "Please insert "<< spacing << " spacing in direction " << axisnames[axis];
            myRead<>(ss.str(), AuxThickness, std::string("Please insert a valid input"));
            ss.str("");
            ss.clear();
            thickness[axis].push_back(AuxThickness);
        };
        this->cell_total=this->cell_total*cell_number[axis];
    };
    defined = true;
};

void Mesh::appear(std::string& _timestamp, int stencil[2]){
    if(_timestamp=="" && this->defined){
        int index = 0;
        for(int axisz=0;axisz<cell_number[2];axisz++){
            for(int axisy=0;axisy<cell_number[1];axisy++){
                for(int axisx=0;axisx<cell_number[0];axisx++){
                    Cell myCell = Cell(index);
                    myCell.setVolume(thickness[2][axisz], thickness[1][axisy], thickness[0][axisx]);
                    myCell.setNumeration3D(0,axisx);
                    myCell.setNumeration3D(1,axisy);
                    myCell.setNumeration3D(2,axisz);
                    cells.push_back(myCell);
                    index++;
                };
            };
        };
        int face_index=0;
        for(auto Celli : cells){
            
            auto LocalNumeration = Celli.getNumeration3D();
            
            for(int axis=0; axis<3; axis++){
                for(int direction=0;direction<2;direction++){
                    LocalNumeration[axis]=LocalNumeration[axis]+stencil[direction];
                    int neighbor_cell=listCell(LocalNumeration);
                    
                    LocalNumeration[axis]=LocalNumeration[axis]-stencil[direction];
                    
                    if(neighbor_cell > -1){
                        
                        Celli.Number_of_active_faces=Celli.Number_of_active_faces+1;
                        Face face = Face();
                        switch(axis){
                        case 0:
                            face.setArea(thickness[1][LocalNumeration[1]]*thickness[2][LocalNumeration[2]]);
                            face.setOrientation(0);
                            break;
                        case 1:
                            face.setArea(thickness[0][LocalNumeration[0]]*thickness[2][LocalNumeration[2]]);
                            face.setOrientation(1);
                            break;
                        case 2:
                            face.setArea(thickness[0][LocalNumeration[0]]*thickness[1][LocalNumeration[1]]);
                            face.setOrientation(2);
                            break;
                        }
                        face.setNeighbor(cells[neighbor_cell]);
                        face_index++;
                        face.setIndex(face_index);
                        Celli.pushFace(face);
                    }
                    
                }
            }
        }
        
    };
};

int Mesh::listCell(int posx, int posy, int posz){
    int index;
    if (!(posx < 0 || posy < 0 && posz < 0 || posx >= cell_number[0] || posy >= cell_number[1] || posz >= cell_number[2])){
        index = posx + posy*cell_number[0] + posz*cell_number[0]*cell_number[1];
	return index;
    }else{
        return -1;
    }
}

int Mesh::listCell(std::vector<int> _Numeration){
    int index;
    if (!(_Numeration[0] < 0 || _Numeration[1] < 0 || _Numeration[2] > 0 || _Numeration[0] >= cell_number[0] || _Numeration[1] >= cell_number[1] || _Numeration[2] >= cell_number[2])){
        index = _Numeration[0] + _Numeration[1]*cell_number[0] + _Numeration[2]*cell_number[0]*cell_number[1];
	return index;
    }else{
        return -1;
    }
    
}

#endif /* MESH_H */
