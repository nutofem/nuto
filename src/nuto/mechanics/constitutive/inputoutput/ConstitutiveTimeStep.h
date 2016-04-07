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
class ConstitutiveTimeStep: public ConstitutiveVector<Eigen::Dynamic>
{
public:
    ConstitutiveTimeStep(int rNumTimeSteps)
    {
        this->resize(rNumTimeSteps);
        this->SetZero();
    }


    //! @brief sets the current time step and shifts the existing time steps by one
    //! @param rCurrentTimeStep ... new current time step
    void SetCurrentTimeStep(double rCurrentTimeStep)
    {
        for (int i = this->rows()-1; i > 0; --i)
            (*this)[i] = (*this)[i-1];
        (*this)[0] = rCurrentTimeStep;
    }

    int GetNumTimeSteps() const
    {
        return this->rows();
    }

};

} /* namespace NuTo */
