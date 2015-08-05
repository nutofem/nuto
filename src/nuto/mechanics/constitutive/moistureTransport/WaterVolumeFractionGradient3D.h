#ifndef WATERVOLUMEFRACTIONGRADIENT3D_H
#define WATERVOLUMEFRACTIONGRADIENT3D_H


#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... water volume fraction gradient 2D
//! @author Volker Hirthammer, BAM
//! @date July 2015
class WaterVolumeFractionGradient3D : public ConstitutiveInputBase, public FullVector<double,3>
{
public:
    WaterVolumeFractionGradient3D();

    virtual const WaterVolumeFractionGradient3D& GetWaterVolumeFractionGradient3D() const override;
};

} // Namespace NuTo

#endif // WATERVOLUMEFRACTIONGRADIENT3D_H
