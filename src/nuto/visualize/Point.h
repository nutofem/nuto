// $Id$

#ifndef POINT_H_
#define POINT_H_
#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/visualize/VisualizeDataBase.h"

namespace NuTo
{
//! @brief ... point for visualization
//! @author Stefan Eckardt, ISM
//! @date November 2009
class Point
{
public:
    //! @brief ... constructor
    //! @param rCoordinates ... point coordinates
    Point(const double* rCoordinates);

    //! @brief ... add scalar data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    void AddDataScalar(unsigned int rDataIndex);

    //! @brief ... add vector (length 3) data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    void AddDataVector(unsigned int rDataIndex);

    //! @brief ... set scalar data
    //! @param rDataIndex ... data index
    //! @param rData ... scalar data
    void SetDataScalar(unsigned int rDataIndex, double rData);

    //! @brief ... set vector data
    //! @param rDataIndex ... data index
    //! @param rData ... vector data
    void SetDataVector(unsigned int rDataIndex, double rData[3]);

    //! @brief ... add tensor (3x3) data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    void AddDataTensor(unsigned int rDataIndex);

    //! @brief ... add field data
    //! @param rDataIndex ... index in data vector (zero based indexing)
    //! @param rNumData ... number of data
    void AddDataField(unsigned int rDataIndex, unsigned int rNumData = 0);

    //! @brief ... get point coordinates
    //! @return ... pointer to point coordinates
    inline const double* GetCoordinates() const
    {
        return this->mCoordinates;
    }

    const VisualizeDataBase* GetData(unsigned int rDataIndex) const;
private:
    //! @brief ... point coordinates
    double mCoordinates[3];

    //! @brief ... vector of point data
    boost::ptr_vector<VisualizeDataBase> mData;
};

}

#endif // POINT_H_ 
