#ifndef VTKMESH_H
#define VTKMESH_H

#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>


class VTKMesh{
 private:

    std::string _filename;

    std::ostringstream _pvd;

    const std::string _header = "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n    <Collection>\n";

    const std::string _ending = "    </Collection>\n</VTKFile>";
    
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();

    vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    
 public:

    VTKMesh(){};

    void filename(const char* filename){_filename = filename;};
    
    template <typename Mesh> void set(const Mesh& mesh){
        auto cell_number = mesh.cellNumber();

        double x=0;
        double y=0;
        double z=0;
        
        for(int k=0; k<=cell_number[2]; ++k){
            
            y=0;
            
            for(int j=0; j<=cell_number[1]; ++j){
                
                x=0;
                
                for(int i=0; i<=cell_number[0]; ++i){

                    points->InsertNextPoint(x, y, z + mesh.top(j,i));
                    x += mesh.thickness(0,i);
                    
                };
                
                y += mesh.thickness(1,j);

            };
            z += mesh.thickness(2,k);
        };

        points->Print(std::cout);

        structuredGrid->SetDimensions(cell_number[0]+1,cell_number[1]+1,cell_number[2]+1);
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

        std::ostringstream ss;

        std::ofstream pvdfile;

        pvdfile.open("./VTKResult/structuredGrid.pvd",ios::trunc);

        ss << "structuredGrid" << "_" << timestep << ".vts";
        _pvd << "    <DataSet timestep=\""<<timestep<<"\" group=\"\" part=\"0\" file=\""<< ss.str() <<"\"/>\n";

        pvdfile << _header <<_pvd.str() << _ending;

        pvdfile.close();
        //writer->SetDataModeToAscii();
        writer->SetFileName(("./VTKResult/"+ss.str()).c_str());
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(structuredGrid);
#else
        writer->SetInputData(structuredGrid);
#endif
        writer->Write();
    };
};

#endif /*VTKMESH_H*/
