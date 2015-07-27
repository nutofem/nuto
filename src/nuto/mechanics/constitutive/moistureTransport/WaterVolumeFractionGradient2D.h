#ifndef WATERVOLUMEFRACTIONGRADIENT2D_H
#define WATERVOLUMEFRACTIONGRADIENT2D_H


#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... water volume fraction gradient 2D
//! @author Volker Hirthammer, BAM
//! @date July 2015
class WaterVolumeFractionGradient2D : public ConstitutiveInputBase, public FullVector<double,2>
{
public:
    WaterVolumeFractionGradient2D();

    virtual const WaterVolumeFractionGradient2D& GetWaterVolumeFractionGradient2D() const override;
};

} // Namespace NuTo

#endif // WATERVOLUMEFRACTIONGRADIENT2D_H
