/*
 * ImplEx.cpp
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#include "mechanics/timeIntegration/ImplEx.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputDummy.h"
#include <algorithm>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

NuTo::ImplEx::ImplEx(StructureBase* rStructure) : ImplicitExplicitBase(rStructure)
{
    mExtrapolationErrorThreshold = -1.;
}

NuTo::ImplEx::~ImplEx()
{}

void NuTo::ImplEx::ExtrapolateStaticData(const ConstitutiveTimeStep& rTimeStep)
{
    // Setup input list with the timeStep and the option to use EULER_FORWARD
    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
            eCalculateStaticData::EULER_FORWARD);
    input[Constitutive::eInput::TIME_STEP] = std::make_unique<ConstitutiveTimeStep>(rTimeStep);


    // Setup output map to store the static data
    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalUpdateStaticData;
    StructureOutputDummy dummy;
    evalUpdateStaticData[eStructureOutput::UPDATE_STATIC_DATA] = &dummy;

    mStructure->Evaluate(input, evalUpdateStaticData);
}

bool NuTo::ImplEx::CheckExtrapolationAndAdjustTimeStep()
{
    if (mExtrapolationErrorThreshold == -1.)
        throw MechanicsException(__PRETTY_FUNCTION__, "Define the extrapolation threshold first.");

    auto errors = mStructure->ElementTotalGetStaticDataExtrapolationError();

    std::partial_sort(errors.begin(), errors.begin() + 5, errors.end(), std::greater<double>());
    double extrapolationError = errors[4];
    std::cout << errors[0] << '\t' << errors[1] << '\t' << errors[2] << std::endl;


    double newTimeStep = std::sqrt(mExtrapolationErrorThreshold / extrapolationError ) * mTimeStep;

    newTimeStep = std::min(newTimeStep, mTimeStep * 1.3);
    newTimeStep = std::max(newTimeStep, mTimeStep / 1.3);
    mTimeStep = newTimeStep;


//    std::cout << "relative Extrapolation error (>1 is bad): " << extrapolationError / mExtrapolationErrorThreshold << std::endl;
//    bool extrapolationIsOK = extrapolationError < mExtrapolationErrorThreshold;
//    if (not extrapolationIsOK && not mForceAcceptOfNextSolution)
//    {
//        mForceAcceptOfNextSolution = true;
//        return false;
//    }
//
//    // accept solution, either the extrapolation is OK or this solution is forced to be accepted (somehow continue the time stepping)
//
//    mForceAcceptOfNextSolution = false;
    return true;
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ImplEx::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ImplEx::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of ImplEx" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase);

    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of ImplEx" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION
