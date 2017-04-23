// $Id$
#pragma once

#include <string>
#include <vector>
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
    int AddCell(CellBase cell);
    
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
    std::vector<Point> mPoints;

    //! @brief ... vector of cells
    std::vector<CellBase> mCells;

    //! @brief ... vector of point data field names
    std::vector<std::string> mPointDataNames;

    //! @brief ... vector of cell data field names 
    std::vector<std::string> mCellDataNames;

    //! @brief ... check if point ids are defined, throws if not
    //! @param pointIds ... point ids
    void CheckPoints(std::vector<int> pointIds) const;

    //! @brief ... get point data index from identifier, throws if not found
    //! @param name ... identifier
    int GetPointDataIndex(const std::string& name) const;

    //! @brief ... get point data index from identifier, throws if not found
    //! @param name ... identifier
    int GetCellDataIndex(const std::string& name) const;
};

}

