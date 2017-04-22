// $Id$
#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include <string>
#include <map>
#include <eigen3/Eigen/Core>

namespace NuTo
{

class CellBase;
class Point;
class VisualizeDataType;

//! @brief ... visualization of unstructured grids
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeUnstructuredGrid
{
public:

    //! @brief ... export to Vtk datafile
    //! @param rFilename ... filename
    void ExportVtkDataFile(const std::string& rFilename) const;

    //! @brief ... export to Vtu datafile
    //! @param rFilename ... filename
    void ExportVtuDataFile(const std::string& rFilename) const;

    //! @brief ... add Point to unstructured grid
    //! @param rCoordinates ... point coordinates
    //! @return ... point id
    int AddPoint(Eigen::Vector3d coordinates);

    //! @brief ... add vertex cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddVertexCell(const unsigned int* rPoints);

    //! @brief ... add line cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddLineCell(const unsigned int* rPoints);

    //! @brief ... add triangle cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddTriangleCell(const unsigned int* rPoints);

    //! @brief ... add quadrilateral cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddQuadCell(const unsigned int* rPoints);

    //! @brief ... add tetraeder cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddTetraCell(const unsigned int* rPoints);

    //! @brief ... add pyramid cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddPyramidCell(const unsigned int* rPoints);

    //! @brief ... add hexahedron cell
    //! @param rPoints ... point id's (zero based indexing)
    //! @return ... cell identifier (zero based indexing)
    unsigned int AddHexahedronCell(const unsigned int* rPoints);

    //! @brief ... define point data 
    //! @param name ... name of the data field
    void DefinePointData(std::string name);

    //! @brief ... define cell data 
    //! @param name ... name of the data field
    void DefineCellData(std::string name);

    //! @brief ... set scalar point data
    //! @param pointIndex ... point index
    //! @param name ... name of the data field
    //! @param data ... scalar data
    void SetPointData(int pointIndex, const std::string& name, double data);

    //! @brief ... set point data
    //! @param pointIndex ... point index
    //! @param name ... name of the data field
    //! @param data ... data
    void SetPointData(int pointIndex, const std::string& name, Eigen::VectorXd data); 

    //! @brief ... set scalar cell data
    //! @param cellIndex ... cell index
    //! @param name ... name of the data field
    //! @param data ... scalar data
    void SetCellData(int cellIndex, const std::string& name, double data);

    //! @brief ... set cell data
    //! @param cellIndex ... cell index
    //! @param name ... name of the data field
    //! @param data ... data
    void SetCellData(int cellIndex, const std::string& name, Eigen::VectorXd data); 


private:
    //! @brief ... vector of points
    boost::ptr_vector<Point> mPoints;

    //! @brief ... vector of cells
    boost::ptr_vector<CellBase> mCells;

    //! @brief ... vector of point data field names
    std::vector<std::string> mPointDataNames;

    //! @brief ... vector of cell data field names 
    std::vector<std::string> mCellDataNames;

    //! @brief ... check if points are defined
    //! @param rNumPoints ... number of points
    //! @param rPoints ... point id's
    void CheckPoints(const unsigned int rNumPoints, const unsigned int *rPoints) const;

    //! @brief ... get point data index from identifier
    //! @param name ... identifier
    int GetPointDataIndex(const std::string& name) const;

    //! @brief ... get point data index from identifier
    //! @param name ... identifier
    int GetCellDataIndex(const std::string& name) const;
};

}

