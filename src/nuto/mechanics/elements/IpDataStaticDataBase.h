// $Id$ 
#pragma once

#include "nuto/mechanics/elements/IpDataBase.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace NuTo
{
class ConstitutiveStaticDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataStaticDataBase : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataStaticDataBase();

    virtual ~IpDataStaticDataBase();

    ConstitutiveStaticDataBase* GetStaticData(int rTimeStep = 0) override;

    const ConstitutiveStaticDataBase* GetStaticData(int rTimeStep = 0) const override;

    void SetStaticData(ConstitutiveStaticDataBase* rStaticData) override;

    IpDataStaticDataBase& GetStaticDataBase()
    {
        return *this;
    }

    const IpDataStaticDataBase& GetStaticDataBase() const
    {
        return *this;
    }

    //! @brief copies the "current" ([0]) static data rNumAdditionalStaticData times
    void AllocateAdditionalStaticData(int rNumAdditionalStaticData) override;

    void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive) override;

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName, double rMacroLength, double rCoordinates[2], std::string rIPName);

    //! @brief sets the fine scale parameter
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(const std::string& rName, double rParameter);

    //! @brief sets the fine scale parameter
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(const std::string& rName, std::string rParameter);

    //! @brief puts current static data to previous static data, previous to pre-previous, etc.
    //! the current data is copied to the previous data, all others are moved
    void SaveStaticData();

    //! @brief puts previous static data to current static data, pre-previous to previous, etc.
    void RestoreStaticData();

    //! @brief returns the total number of allocated ip data sets
    int GetNumStaticData() const;

#ifdef ENABLE_VISUALIZE
	//Visualize for all integration points the fine scale structure
	void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
			const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;
#endif

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:

    boost::ptr_vector<ConstitutiveStaticDataBase> mStaticData;
};
}
