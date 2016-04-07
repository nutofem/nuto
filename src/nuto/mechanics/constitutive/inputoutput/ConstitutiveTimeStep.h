/*
 * ConstitutiveTimeStep.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"

namespace NuTo
{

//! @brief class to store the time step dT.
//! dT[0] = current time step
//! dT[1] = previous time step
//! dT[2] = pre-previous time step
//! ...
template <int TNum>
class ConstitutiveTimeStep: public ConstitutiveVector<TNum>
{
public:
    ConstitutiveTimeStep()
    {
        this->SetZero();
    }


    //! @brief sets the current time step and shifts the existing time steps by one
    //! @param rCurrentTimeStep ... new current time step
    void SetCurrentTimeStep(double rCurrentTimeStep)
    {
        for (int i = TNum-1; i > 0; --i)
            (*this)[i] = (*this)[i-1];
        (*this)[0] = rCurrentTimeStep;
    }

};

} /* namespace NuTo */
