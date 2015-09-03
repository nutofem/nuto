//============================================================================
// Name        : InterfaceSlip.cpp
// Author      : Philip Huschke
// Version     : 26 Aug 2015
// Copyright   :
// Description : Constitutive input for the interface element proposed by Goodman et al.
//============================================================================

#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{

class InterfaceSlip: public ConstitutiveInputBase
{

public:
    //! @brief Constructor
    InterfaceSlip();

    //! @brief  Get interface slip
    //! @return Pointer to the object
    const NuTo::InterfaceSlip& GetInterfaceSlip() const override;

    const Eigen::VectorXd GetInterfaceSlipVector() const;

    //! @brief Get interface slip
    //! @param  rInterfaceSlip: interface slip
    void GetInterfaceSlip(Eigen::VectorXd& rInterfaceSlip) const;

    //! @brief Set interface slip
    //! @param  rInterfaceSlip: interface slip
    void SetInterfaceSlip(const Eigen::VectorXd& rInterfaceSlip);

private:
    Eigen::VectorXd mInterfaceSlip;
};

}

