#include "mechanics/constitutive/mechanics/InterfaceSlip.h"

NuTo::InterfaceSlip::InterfaceSlip()
{
    mInterfaceSlip.setZero();
}

const NuTo::InterfaceSlip& NuTo::InterfaceSlip::GetInterfaceSlip() const
{
    return *this;
}

const Eigen::VectorXd NuTo::InterfaceSlip::GetInterfaceSlipVector() const
{
    return mInterfaceSlip;
}

void NuTo::InterfaceSlip::GetInterfaceSlip(Eigen::VectorXd& rInterfaceSlip) const
{
    rInterfaceSlip = GetInterfaceSlipVector();
}


void NuTo::InterfaceSlip::SetInterfaceSlip(const Eigen::VectorXd& rInterfaceSlip)
{
    mInterfaceSlip = rInterfaceSlip;
}

