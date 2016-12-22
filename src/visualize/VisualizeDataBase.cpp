// $Id$

#include "visualize/VisualizeDataBase.h"

namespace NuTo
{
std::ostream& operator<<(std::ostream& os, const VisualizeDataBase& rData)
{
    return rData.Output(os);
}
}
