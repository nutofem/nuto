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
#include "base/serializeStream/SerializeStreamOut.h"
#include "mechanics/timeIntegration/ImplExCallback.h"

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
    mImplExCallback = new ImplExCallback();
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

    double extrapolationError = mStructure->ElementTotalGetStaticDataExtrapolationError();

    mTimeStep = mImplExCallback->GetNewTimeStep(extrapolationError, mExtrapolationErrorThreshold, mTimeStep);
    return mImplExCallback->AcceptSolution(extrapolationError, mExtrapolationErrorThreshold);
}

void NuTo::ImplEx::SetImplExCallback(NuTo::ImplExCallback* r)
{
    if (mImplExCallback != nullptr)
        delete mImplExCallback;
    mImplExCallback = r;
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
