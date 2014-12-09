#ifndef WATERPHASEFRACTION_H
#define WATERPHASEFRACTION_H

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... water phase fraction
//! @author Volker Hirthammer, BAM
//! @date December 2014
class WaterPhaseFraction : public ConstitutiveInputBase, public FullVector<double,1>
{
public:
    WaterPhaseFraction();

    virtual const WaterPhaseFraction& GetWaterPhaseFraction() const override;

};

} // Namespace NuTo



#endif // WATERPHASEFRACTION_H
