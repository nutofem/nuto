#pragma once


#include "mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"

namespace NuTo
{

class LinearElasticInhomogeneous: public LinearElasticEngineeringStress
{
public:

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<LinearElasticInhomogeneous>>(*this);
    }

    void UpdateParameters(Eigen::VectorXd coordinates);

    void SetYoungsModulus(std::function<double(Eigen::VectorXd)> E);

    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const override;

    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                          const ConstitutiveOutputMap& rConstitutiveOutput);

protected:
    std::function<double(Eigen::VectorXd)> mE;

};

}
