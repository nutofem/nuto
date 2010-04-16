#include <string>
#include "nuto/metamodel/MetamodelException.h"

namespace NuTo
{
//! @brief ... return error message of the exception
//! @return ... error message
std::string MetamodelException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module MetaModel\n"+message);

    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
MetamodelException* MetamodelException::Clone()
{
    return new MetamodelException(*this);
}
} //namespace NuTo
