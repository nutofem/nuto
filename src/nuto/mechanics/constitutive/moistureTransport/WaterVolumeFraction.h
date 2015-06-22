#ifndef WATERVOLUMEFRACTION_H
#define WATERVOLUMEFRACTION_H

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... water volume fraction
//! @author Volker Hirthammer, BAM
//! @date December 2014
class WaterVolumeFraction : public ConstitutiveInputBase, public FullVector<double,1>
{
public:
    WaterVolumeFraction();

    virtual const WaterVolumeFraction& GetWaterVolumeFraction() const override;

};

} // Namespace NuTo



#endif // WATERVOLUMEFRACTION_H
