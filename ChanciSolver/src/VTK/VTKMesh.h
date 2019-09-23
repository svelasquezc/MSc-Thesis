#ifndef VTKMESH_H
#define VTKMESH_H

#include <string>
#include <sstream>
#include <vector>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>


class VTKMesh{
 private:

    std::string filename;
    
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();

    vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    
 public:

    VTKMesh(){};
    
    template <typename Mesh> void set(const Mesh& mesh){
        auto cell_number = mesh.cellNumber();

        double x=0;
        double y=0;
        double z=0;
        
        for(int k=0; k<=cell_number[2]; ++k){

            for(int j=0; j<=cell_number[1]; ++j){

                for(int i=0; i<=cell_number[0]; ++i){

                    points->InsertNextPoint(x, y, z + mesh.top(j,i));
                    x += mesh.thickness(0,i);
                    
                };
                
                y += mesh.thickness(1,j);

            };
            z += mesh.thickness(2,k);
        };

        structuredGrid->SetDimensions(cell_number[0],cell_number[1],cell_number[2]);
        structuredGrid->SetPoints(points);
        
    };

    void appendScalar(const std::string name, const std::vector<double> all_values){
        auto scalar = vtkSmartPointer<vtkDoubleArray>::New();
        scalar->SetName(name.c_str());
        scalar->SetNumberOfComponents(1);
        for(vtkIdType i = 0; i< all_values.size(); ++i){
            scalar->InsertNextValue(all_values[i]);
        }
        structuredGrid->GetCellData()->AddArray(scalar);
        structuredGrid->Modified();
    };
    
    void write(double timestep){

        ostringstream ss; 
        
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

        ss << filename << "_" << timestep << ".vts";
        
        writer->SetFileName(ss.str().c_str());
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(structuredGrid);
#else
        writer->SetInputData(structuredGrid);
#endif
        writer->Write();
    };
};

#endif /*VTKMESH_H*/
