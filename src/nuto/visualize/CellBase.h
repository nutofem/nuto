// $Id$

#ifndef CELLBASE_H_
#define CELLBASE_H_
#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/visualize/VisualizeDataBase.h"
#include "nuto/visualize/VisualizeDataType.h"

namespace NuTo
{
//! @brief ... base class for visualization cells
//! @author Stefan Eckardt, ISM
//! @date November 2009
class CellBase
{
public:
    enum eCellTypes
    {
        LINE,
        TRIANGLE,
        QUAD,
        TETRAEDER,
        HEXAHEDRON,
        PYRAMID
    };

    //! @brief constructor
    //! @param rDataTypes ... data type definitions
    CellBase(const std::vector<VisualizeDataType>& rDataTypes);

    //! @brief ... return number of cell points
    //! @return ... number of cell points
    virtual unsigned int GetNumPoints() const = 0;

    //! @brief ... return point id's
    //! @return ... array of point id's
    virtual const unsigned int* GetPoints() const = 0;

    //! @brief ... returns the corresponding Vtk cell type
    //! @return ... Vtk cell type
    virtual unsigned int GetVtkCellType() const = 0;

    //! @brief ... add scalar data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    void AddDataScalar(unsigned int rDataIndex);

    //! @brief ... add vector (length 3) data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    void AddDataVector(unsigned int rDataIndex);

    //! @brief ... add tensor (3x3) data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    void AddDataTensor(unsigned int rDataIndex);

    //! @brief ... set tensor data
    //! @param rDataIndex ... data index
    //! @param rData ... tensor data
    void SetDataTensor(unsigned int rDataIndex, double rData[9]);

    //! @brief ... set tensor data
    //! @param rDataIndex ... data index
    //! @param rData ... tensor data
    void SetDataScalar(unsigned int rDataIndex, double rData);

    //! @brief ... set vector data
    //! @param rDataIndex ... data index
    //! @param rData ... tensor data
    void SetDataVector (unsigned int rDataIndex, double rData[3]);

    //! @brief ... add field data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    //! @param rNumData ... number of data
    void AddDataField(unsigned int rDataIndex, unsigned int rNumData = 0);

    const VisualizeDataBase* GetData(unsigned int rDataIndex) const;
protected:
    boost::ptr_vector<VisualizeDataBase> mData;
};

}

#endif // CELLBASE_H_ 
