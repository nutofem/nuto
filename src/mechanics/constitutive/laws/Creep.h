#pragma once


#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/staticData/DataCreep.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"


namespace NuTo
{
class Creep : public ConstitutiveBase
{
public:
    typedef Constitutive::StaticData::DataCreep StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<StaticDataType>;

    // Ctor / Dtor
    // -----------
    Creep();
    Creep(const Creep&) = default;
    Creep(Creep&&) = default;
    ~Creep() = default;

    // assignment operator
    // -------------------
    Creep& operator=(const Creep&) = delete;
    Creep& operator=(Creep&&) = delete;


    // evaluate

    //! @brief ... evaluate the constitutive relation
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    //! @param rStaticData ... static data
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput,
                  Data& rStaticData);


    // pure virtual overrides

    virtual std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<Creep>>(*this, StaticDataType());
    }

    //! @brief Determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput Desired constitutive outputs
    //! @return constitutive inputs needed for the evaluation
    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const override;

    //! @brief Determines which submatrices of a multi-doftype problem can be solved by the constitutive law.
    //! @param rDofRow Row DOF.
    //! @param rDofCol Column DOF.
    //! @param rTimeDerivative Time derivative.
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or
    //! stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const override;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrown
    virtual void CheckParameters() const override;


    // setter
    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter identifier, double value) override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param identifier ... Enum to identify the requested parameter
    //! @param value ... new value for requested variable
    virtual void SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter identifier,
                                              Eigen::VectorXd value) override;


private:
    //! @brief Calculates the stiffness matrix
    //! @return Stiffness matrix
    template <int TDim>
    Eigen::MatrixXd ExponentialAlgorithmCalculateStiffnessMatrix(double delta_t,
                                                                 const ConstitutiveInputMap& rConstitutiveInput) const;

    //! @brief ... Calculates the beta parameter as described in TODO: INSERT REFERENCE
    double ExponentialAlgorithmCalculateBeta(unsigned int index, double delta_t) const;

    //! @brief ... Calculates the beta parameter as described in TODO: INSERT REFERENCE
    double ExponentialAlgorithmCalculateLambda(unsigned int index, double delta_t) const;

    //! @brief ... Calculates the beta parameter as described in TODO: INSERT REFERENCE
    double ExponentialAlgorithmCalculateChainStiffness(double delta_t) const;


protected:
    Eigen::VectorXd mKC_E = {};
    Eigen::VectorXd mKC_D = {};
    Eigen::VectorXd mKC_T = {};
    double mE = 0.0;
    double mNu = 0.0;
};


} // namespace NuTo
