// $Id$

#pragma once

#include <iostream>

#include "visualize/VisualizeDataBase.h"

namespace NuTo
{

//! @brief ... field (variable size) data for visualization
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeDataField: public NuTo::VisualizeDataBase
{
public:
    // constructor
    VisualizeDataField(unsigned int rNumData);

    // destructor
    ~VisualizeDataField();

    //! @brief ... get data type
    //! @return ... visualize data type
    NuTo::eVisualizeDataType GetDataType() const;

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
    unsigned int mNumData;
    double* mData;

    //! @brief ... create output stream
    //! @param os ... output stream
    //! @return ... output stream
    std::ostream& Output(std::ostream& os) const;
};
}


