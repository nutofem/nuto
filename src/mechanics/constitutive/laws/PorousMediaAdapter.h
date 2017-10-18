#pragma once

#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/PDEs/PorousMedium.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{

class PorousMediaAdapter : public ConstitutiveBase
{
public:
    PorousMediaAdapter(double a, double b, double c, double d)
        : mProperLaw(a, b, c, d)
    {
    }

    const PorousMedium& GetProperLaw() const
    {
        return mProperLaw;
    }

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<PorousMediaAdapter>>(*this);
    }

    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const override
    {
        return ConstitutiveInputMap();
    }

    bool CheckDofCombinationComputable(Node::eDof dofRow, Node::eDof dofCol, int timeDerivative) const override
    {
        if (timeDerivative == 2)
            return false;
        else if (dofRow == Node::eDof::CAPILLARY_PRESSURE and dofCol == Node::eDof::CAPILLARY_PRESSURE)
            return true;
        else if (dofRow == Node::eDof::CAPILLARY_PRESSURE and dofCol == Node::eDof::GAS_PRESSURE)
            return true;
        else if (dofRow == Node::eDof::GAS_PRESSURE and dofCol == Node::eDof::CAPILLARY_PRESSURE)
            return true;
        else if (dofRow == Node::eDof::GAS_PRESSURE and dofCol == Node::eDof::GAS_PRESSURE)
            return true;
        else
            return false;
    }

    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
    {
    }

    Constitutive::eConstitutiveType GetType() const override
    {
        return Constitutive::eConstitutiveType::WAHWAHWAH;
    }

    bool HaveTmpStaticData() const override
    {
        return false;
    }

    void CheckParameters() const override {}

private:
    PorousMedium mProperLaw;
};
}
