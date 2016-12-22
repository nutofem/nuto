#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include "mechanics/constitutive/staticData/DataAdditiveInputImplicit.h"
#include "mechanics/constitutive/staticData/DataContainer.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"
#include <type_traits>

namespace NuTo
{
class AdditiveInputImplicit;
namespace Constitutive
{

class IPAdditiveInputImplicit: public IPConstitutiveLawBase
{
public:

    typedef StaticData::DataAdditiveInputImplicit StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<StaticDataType>;

    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPAdditiveInputImplicit(AdditiveInputImplicit& rLaw, const Data &rData);

    IPAdditiveInputImplicit(const IPAdditiveInputImplicit& rOther);

    std::unique_ptr<IPConstitutiveLawBase> Clone() const override;

    ConstitutiveBase& GetConstitutiveLaw() const override;

    //! @brief allocates rNum additional static data
    //! @param rNum number of addtional static data
    void AllocateAdditional(int rNum) override;

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    void ShiftToPast() override;

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    void ShiftToFuture() override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) override;

protected:
    template<int TDim>
    eError AdditiveInputImplicitEvaluate(const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput);


    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return AdditiveInputImplicitEvaluate<1>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return AdditiveInputImplicitEvaluate<2>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return AdditiveInputImplicitEvaluate<3>(rConstitutiveInput, rConstitutiveOutput);
    }

private:

    template<int TDim>
    void CalculateGlobalOutputs(const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                std::vector<NuTo::ConstitutiveInputMap> &rLocalInputMapVec,
                                std::vector<ConstitutiveOutputMap> &rLocalOutputMapVec);

    template<int TDim>
    void CreateLocalInAndOutputMaps(const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                    const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                    std::vector<NuTo::ConstitutiveInputMap> &rLocalInputMapVec,
                                    std::vector<ConstitutiveOutputMap> &rLocalOutputMapVec);

    AdditiveInputImplicit& mLaw;
    Data mData;
    boost::ptr_vector<IPConstitutiveLawBase> mSublawIPs;
};


} // namespace Constitutive
} // namespace NuTo
