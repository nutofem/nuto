// $Id$

#pragma once

#include <iostream>

namespace NuTo
{


enum class eVisualizeDataType;

//! @brief ...
//! @author Stefan Eckardt, ISM
//! @date 23.11.2009
class VisualizeDataBase
{
public:
    virtual ~VisualizeDataBase(){};

    //! @brief ... get data type
    //! @return ... visualize data type
    virtual NuTo::eVisualizeDataType GetDataType() const = 0;

    //! @brief ... get number of data
    //! @return ... number of data
    virtual unsigned int GetNumData() const = 0;

    //! @brief ... get data
    //! @return ... pointer to data array
    virtual const double* GetData() const = 0;

    //! @brief ... set data
    //! @param rData ... data array
    virtual void SetData(const double* rData) = 0;

    //! @brief ... stream operator
    //! @param os ... output stream
    //! @param rData ... data object
    friend std::ostream& operator<<(std::ostream& os, const VisualizeDataBase& rData);

protected:
    //! @brief ... create output stream
    //! @param os ... output stream
    //! @return ... output stream
    virtual std::ostream& Output(std::ostream& os) const = 0;
};

}

