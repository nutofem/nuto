#ifndef RELATIVEHUMIDITYGRADIENT2D_H
#define RELATIVEHUMIDITYGRADIENT2D_H


#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... water volume fraction
//! @author Volker Hirthammer, BAM
//! @date December 2014
class RelativeHumidityGradient2D : public ConstitutiveInputBase, public FullVector<double,2>
{
public:
    RelativeHumidityGradient2D();

    virtual const RelativeHumidityGradient2D& GetRelativeHumidityGradient2D() const override;

};

} // Namespace NuTo

#endif // RELATIVEHUMIDITYGRADIENT2D_H
