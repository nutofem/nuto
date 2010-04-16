// $Id$

#ifndef VISUALIZEDATABASE_H_
#define VISUALIZEDATABASE_H_
#include <iostream>
#include "nuto/visualize/VisualizeDataType.h"
namespace NuTo
{
//! @brief ...
//! @author Stefan Eckardt, ISM
//! @date 23.11.2009
class VisualizeDataBase
{
public:
    //! @brief ... get data type
    //! @return ... visualize data type
    virtual NuTo::VisualizeDataType::eDataType GetDataType() const = 0;

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

private:
    //! @brief ... create output stream
    //! @param os ... output stream
    //! @return ... output stream
    virtual std::ostream& Output(std::ostream& os) const = 0;
};

}

#endif // VISUALIZEDATABASE_H_
