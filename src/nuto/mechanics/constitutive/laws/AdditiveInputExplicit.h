#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataMultipleConstitutiveLaws.h"

#include <set>
#include <vector>

namespace NuTo
{

class AdditiveInputExplicit : public ConstitutiveBase
{
public:

    //! @brief constructor
    AdditiveInputExplicit()
        : ConstitutiveBase()
    {
        mComputableDofCombinations.resize(2); //VHIRTHAMTODO ---> Get number time derivatives during construction (as parameter)
    }


    //! @brief ... adds a constitutive law to a model that combines multiple constitutive laws (additive, parallel)
    //! @param rConstitutiveLaw ... additional constitutive law
    //! @param rModiesInput ... enum which defines wich input is modified by a constitutive law.
    virtual void  AddConstitutiveLaw(NuTo::ConstitutiveBase* rConstitutiveLaw, Constitutive::Input::eInput rModiesInput = Constitutive::Input::NONE) override;


    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData1D(const ElementBase* rElement) const override
    {
        return AllocateStaticDataAdditiveInputExplicit<1>(rElement);
    }

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData2D(const ElementBase* rElement) const override
    {
        return AllocateStaticDataAdditiveInputExplicit<2>(rElement);
    }

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData3D(const ElementBase* rElement) const override
    {
        return AllocateStaticDataAdditiveInputExplicit<3>(rElement);
    }



    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow,
                                               Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility( Element::eElementType rElementType) const override
    {
        for(unsigned int i=0; i<mSublaws.size(); ++i)
        {
            if(!mSublaws[i].first->CheckElementCompatibility(rElementType))
                return false;
        }
        assert(mMainLaw!=nullptr);
        return mMainLaw->CheckElementCompatibility(rElementType);
    }

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void CheckParameters() const override
    {
        for(unsigned int i=0; i<mSublaws.size(); ++i)
        {
            mSublaws[i].first->CheckParameters();
        }
        assert(mMainLaw!=nullptr);
        mMainLaw->CheckParameters();
    }




    //! @brief ... evaluate the constitutive relation of every attached constitutive law in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate1D(ElementBase* rElement,
                                           int rIp,
                                           const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return EvaluateAdditiveInputExplicit<1>(rElement,
                                                rIp,
                                                rConstitutiveInput,
                                                rConstitutiveOutput);
    }

    //! @brief ... evaluate the constitutive relation of every attached constitutive law in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate2D(ElementBase* rElement,
                                           int rIp,
                                           const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return EvaluateAdditiveInputExplicit<2>(rElement,
                                                rIp,
                                                rConstitutiveInput,
                                                rConstitutiveOutput);

    }

    //! @brief ... evaluate the constitutive relation of every attached constitutive law relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate3D(ElementBase* rElement,
                                           int rIp,
                                           const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return EvaluateAdditiveInputExplicit<3>(rElement,
                                                rIp,
                                                rConstitutiveInput,
                                                rConstitutiveOutput);
    }



    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    virtual ConstitutiveInputMap GetConstitutiveInputs( const ConstitutiveOutputMap& rConstitutiveOutput,
                                                        const InterpolationType& rInterpolationType) const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override
    {
        return NuTo::Constitutive::ADDITIVE_INPUT_EXPLICIT;
    }

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const override
    {
        for(unsigned int i=0; i<mSublaws.size(); ++i)
        {
            if(mSublaws[i].first->HaveTmpStaticData())
                return true;
        }
        if(mMainLaw!=nullptr)
        {
            return mMainLaw->HaveTmpStaticData();
        }
        return false;
    }

private:

    //! @brief Adds all calculable dof combinations of an attached constitutive law to an internal storage
    //! @param rConstitutiveLaw ... constitutive law
    void AddCalculableDofCombinations(NuTo::ConstitutiveBase* rConstitutiveLaw);

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    template <int TDim>
    ConstitutiveStaticDataBase* AllocateStaticDataAdditiveInputExplicit(const ElementBase* rElement) const;

    //! @brief Applies sublaw dependent modifications to the main laws inputs and the global outputs
    //! @param rMainLawInput: Input map of the main law
    //! @param rConstitutiveOutput: Global output map of this constitutive law
    //! @param rSublawOutput: Output map of a sublaw
    template <int TDim>
    void ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput,
                            const ConstitutiveOutputMap& rConstitutiveOutput,
                            const ConstitutiveOutputMap& rSublawOutput);

    //! @brief Calculates the derivatives that depend on main law and sublaw outputs
    //! @param rConstitutiveOutput: Global output map of this constitutive law
    //! @param rSublawOutputVec: Vector with all the sublaws output maps
    template <int TDim>
    void CalculateDerivatives(const ConstitutiveOutputMap &rConstitutiveOutput,
                              std::vector<ConstitutiveOutputMap> &rSublawOutputVec);

    //! @brief ... evaluate the constitutive relation
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    template <int TDim>
    NuTo::Error::eError EvaluateAdditiveInputExplicit(  ElementBase* rElement,
                                                        int rIp,
                                                        const ConstitutiveInputMap& rConstitutiveInput,
                                                        const ConstitutiveOutputMap& rConstitutiveOutput);


    //! @brief Gets the enum of the sublaw output that is needed to calculate the specified derivative
    //! @param rParameter: Enum of the parameter whose derivative is needed
    //! @param rMainDerivative: The requested global derivative
    //! @return Sublaw derivative enum
    Constitutive::Output::eOutput GetDerivativeEnumSublaw(Constitutive::Output::eOutput rParameter,
                                                          Constitutive::Output::eOutput rMainDerivative) const;


    //! @brief Gets a modified output map for the sublaws depending on the main laws inputs and the specified modified Input (see Member: mModifiedInputs)
    //! @return Modified output map
    template <int TDim>
    ConstitutiveOutputMap GetSublawOutputMap(const ConstitutiveInputMap& rMainLawInputMap,
                                             const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                             unsigned int rSublawIndex) const;






    //! @brief ... pointer to constitutive law that delivers the output
    NuTo::ConstitutiveBase* mMainLaw = nullptr;

    //! @brief ... vector of pairs which store the pointers to input modifying constitutive laws together with the modified inputs enum.
    std::vector<std::pair<NuTo::ConstitutiveBase*,Constitutive::Input::eInput>> mSublaws;

    std::vector<std::set<std::pair<Node::eDof,Node::eDof>>> mComputableDofCombinations;

    //! @brief debug variable to avoid that a constitutive law can be attached after allocation of static data.
    mutable bool mStaticDataAllocated = false;
};

}
