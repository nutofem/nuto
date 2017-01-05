#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <eigen3/Eigen/Core>

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include "base/ErrorEnum.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawBase.h"

namespace NuTo
{
class InterpolationType;
class ElementBase;
class InterpolationType;
class Logger;
enum class eError;
template<typename IOEnum> class ConstitutiveIOMap;
template <class T, int rows> class FullVector;
template <class T, int rows, int cols> class FullMatrix;

namespace Element
{
    enum class eElementType;
}// namespace Element

namespace Constitutive
{
    enum class eConstitutiveParameter;
    enum class eConstitutiveType;
    enum class eInput;
    enum class eOutput;
}// namespace Constitutive

namespace Node
{
    enum class eDof : unsigned char;
}// namespace Node
using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::eInput>;
using ConstitutiveOutputMap = ConstitutiveIOMap<Constitutive::eOutput>;

//! @brief Base class for the constitutive relationship, e.g. material laws.
class ConstitutiveBase
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutiveBase(): mParametersValid(false){};

    ConstitutiveBase(const ConstitutiveBase& ) = default;
    ConstitutiveBase(      ConstitutiveBase&&) = default;

    ConstitutiveBase& operator =(const ConstitutiveBase& ) = default;
    ConstitutiveBase& operator =(      ConstitutiveBase&&) = default;

    //! @brief ... constructor
    virtual ~ConstitutiveBase() = default;

    virtual std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() = 0;

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       const InterpolationType& rInterpolationType) const = 0;

    //! @brief Determines which submatrices of a multi-doftype problem can be solved by the constitutive law.
    //! @param rDofRow Row DOF.
    //! @param rDofCol Column DOF.
    //! @param rTimeDerivative Time derivative.
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
            int rTimeDerivative) const = 0;

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... checks if the constitutive law has a specific parameter
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... true/false
    virtual bool CheckHaveParameter(Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual bool GetParameterBool(Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterBool(Constitutive::eConstitutiveParameter rIdentifier, bool rValue);

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue);

    virtual void SetParameterFunction(std::function<std::array<double, 2>(double)>);

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual Eigen::VectorXd GetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier, Eigen::VectorXd rValue);


    //! @brief checks parameters, throws if the check failed
    static void CheckParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue);

    //! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
    //! @param rRelativeHumidity ... relative humidity
    //! @return ... equilibrium water volume fraction
    virtual double GetEquilibriumWaterVolumeFraction(double rRelativeHumidity, Eigen::VectorXd rCoeffs) const;

    //! @brief ... checks if a constitutive law has an specific output
    //! @return ... true/false
    virtual bool CheckOutputTypeCompatibility(NuTo::Constitutive::eOutput rOutputEnum) const;

    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const = 0;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const = 0;

    //! @brief ... returns whether the parameters of the constitutive relationship are valid or not
    //! @return ...  <B>true</B> if all parameters of the constitutive relationship are valid and <B>false</B> otherwise
    inline bool AreParametersValid() const
    {
        return this->mParametersValid;
    }

    //! @brief ... check if all parameters are valid and modify parameter validity flag
    //! @sa mParametersValid
    void SetParametersValid();

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel, Logger& rLogger) const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const=0;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrown
    virtual void CheckParameters() const = 0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... flag which is <B>true</B> if all parameters of the constitutive relationship are valid and <B>false</B> otherwise
    bool mParametersValid;

};

}
