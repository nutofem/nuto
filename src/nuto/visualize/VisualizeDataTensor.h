// $Id$

#ifndef VISUALIZEDATATENSOR_H_
#define VISUALIZEDATATENSOR_H_
#include <iostream>

#include "nuto/visualize/VisualizeDataBase.h"

namespace NuTo
{

//! @brief ... tensor (size 3x3) data for visualization
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeDataTensor: public NuTo::VisualizeDataBase
{
public:
    // constructor
    VisualizeDataTensor();

    //! @brief ... get data type
    //! @return ... visualize data type
    NuTo::VisualizeDataType::eDataType GetDataType() const;

    //! @brief ... get number of data
    //! @return ... number of data
    unsigned int GetNumData() const;

    //! @brief ... get data
    //! @return ... pointer to data array
    const double* GetData() const;

    //! @brief ... set data
    //! @param rData ... data array
    void SetData(const double* rData);

private:
    double mData[9];

    //! @brief ... create output stream
    //! @param os ... output stream
    //! @return ... output stream
    std::ostream& Output(std::ostream& os) const;
};
}

#endif // VISUALIZEDATATENSOR_H_ 
