// $Id$

#ifndef VISUALIZEUNSTRUCTUREDGRID_H_
#define VISUALIZEUNSTRUCTUREDGRID_H_
#include <boost/ptr_container/ptr_vector.hpp>
#include <map>
#include <string>

#include "nuto/visualize/VisualizeDataType.h"
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/Point.h"
#include "nuto/visualize/CellBase.h"

namespace NuTo
{
//! @brief ... visualization of unstructured grids
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeUnstructuredGrid: public NuTo::VisualizeBase
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
    //! @return ... point identifier (zero based indexing)
    unsigned int AddPoint(const double* rCoordinates);

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

    //! @brief ... define scalar point data at the structure
    //! @param rIdent ... identifier
    void DefinePointDataScalar(const std::string& rIdent);

    //! @brief ... define vector (length 3) point data at the structure
    //! @param rIdent ... identifier
    void DefinePointDataVector(const std::string& rIdent);

    //! @brief ... set scalar point data
    //! @param rPointIndex ... point index
    //! @param rDataIdent ... data identifier
    //! @param rData ... scalar data
    void SetPointDataScalar(unsigned int rPointIndex, const std::string& rDataIdent, double rData);

    //! @brief ... set vector point data
    //! @param rPointIndex ... point index
    //! @param rDataIdent ... data identifier
    //! @param rData ... vector data
    void SetPointDataVector(unsigned int rPointIndex, const std::string& rDataIdent, double rData[3]);

    //! @brief ... define tensor (3x3) point data at the structure
    //! @param rIdent ... identifier
    void DefinePointDataTensor(const std::string& rIdent);

    //! @brief ... define field point data at the structure
    //! @param rIdent ... identifier
    //! @param rNumData ... number of field data components
    void DefinePointDataField(const std::string& rIdent, unsigned int rNumData = 0);

    //! @brief ... define scalar cell data at the structure
    //! @param rIdent ... identifier
    void DefineCellDataScalar(const std::string& rIdent);

    //! @brief ... define vector (length 3) cell data at the structure
    //! @param rIdent ... identifier
    void DefineCellDataVector(const std::string& rIdent);

    //! @brief ... define tensor (3x3) cell data at the structure
    //! @param rIdent ... identifier
    void DefineCellDataTensor(const std::string& rIdent);

    //! @brief ... set scalar cell data
    //! @param rPointIndex ... cell index
    //! @param rDataIdent ... data identifier
    //! @param rData ... scalar data
    void SetCellDataScalar(unsigned int rCellIndex, const std::string& rDataIdent, double rData);

    //! @brief ... set vector cell data
    //! @param rPointIndex ... cell index
    //! @param rDataIdent ... data identifier
    //! @param rData ... vector data
    void SetCellDataVector(unsigned int rCellIndex, const std::string& rDataIdent, double rData[3]);

    //! @brief ... set tensor cell data
    //! @param rPointIndex ... cell index
    //! @param rDataIdent ... data identifier
    //! @param rData ... tensor data
    void SetCellDataTensor(unsigned int rCellIndex, const std::string& rDataIdent, double rData[9]);

    //! @brief ... define field cell data at the structure
    //! @param rIdent ... identifier
    //! @param rNumData ... number of field data components
    void DefineCellDataField(const std::string& rIdent, unsigned int rNumData = 0);

private:
    //! @brief ... vector of points
    boost::ptr_vector<Point> mPoints;

    //! @brief ... vector of cells
    boost::ptr_vector<CellBase> mCells;

    //! @brief ...vector of point data definitions
    std::vector<VisualizeDataType> mPointData;

    //! @brief ... vector of cell data definitions
    std::vector<VisualizeDataType> mCellData;

    //! @brief ... check if points are defined
    //! @param rNumPoints ... number of points
    //! @param rPoints ... point id's
    void CheckPoints(const unsigned int rNumPoints, const unsigned int *rPoints) const;

    //! @brief ... check identifier for point data
    //! @param rIdent ... identifier
    void CheckPointDataIdent(const std::string& rIdent) const;

    //! @brief ... check identifier for cell data
    //! @param rIdent ... identifier
    void CheckCellDataIdent(const std::string& rIdent) const;

    //! @brief ... check identifier for data
    void CheckDataIdent(const std::string& rIdent) const;

    //! @brief ... get point data index from identifier
    //! @param rIdent ... data identifier
    unsigned int GetPointDataIndex(const std::string& rIdent) const;

    //! @brief ... get point data index from identifier
    //! @param rIdent ... data identifier
    unsigned int GetCellDataIndex(const std::string& rIdent) const;
};

}

#endif // VISUALIZEUNSTRUCTUREDGRID_H_ 
