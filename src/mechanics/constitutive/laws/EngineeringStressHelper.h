#pragma once

#include <tuple>
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"

namespace NuTo
{

// forward declarations
class InterpolationType;
template<typename IOEnum> class ConstitutiveIOMap;
template <int TDim> class EngineeringStrain;
template <int TDim> class EngineeringStress;
namespace Constitutive
{
    enum class eInput;
}
using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::eInput>;



class EngineeringStressHelper
{
public:
    //! @brief calculate coefficients of the PLANE_STRESS 2D material matrix
    //! @param rE ... Young's modulus
    //! @param rNu ... Poisson's ratio
    //! @return tuple <C11, C12, C33>
    static std::tuple<double, double, double> CalculateCoefficients2DPlaneStress(double rE, double rNu);

    //! @brief calculate coefficients of the 3D material matrix
    //! @param rE ... Young's modulus
    //! @param rNu ... Poisson's ratio
    //! @return tuple <C11, C12, C33>
    static std::tuple<double, double, double> CalculateCoefficients3D(double rE, double rNu);



    template <int TDim>
    static NuTo::EngineeringStress<TDim> GetStress(const NuTo::EngineeringStrain<TDim>& rElasticStrain,
                                                   double rE,
                                                   double rNu, ePlaneState rPlaneState = ePlaneState::PLANE_STRESS);
};
} /* namespace NuTo */

