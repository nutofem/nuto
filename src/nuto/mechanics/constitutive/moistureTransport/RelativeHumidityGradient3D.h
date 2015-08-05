#ifndef RELATIVEHUMIDITYGRADIENT3D_H
#define RELATIVEHUMIDITYGRADIENT3D_H


#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... water volume fraction
//! @author Volker Hirthammer, BAM
//! @date December 2014
class RelativeHumidityGradient3D : public ConstitutiveInputBase, public FullVector<double,3>
{
public:
    RelativeHumidityGradient3D();

    virtual const RelativeHumidityGradient3D& GetRelativeHumidityGradient3D() const override;

};

} // Namespace NuTo

#endif // RELATIVEHUMIDITYGRADIENT3D_H
