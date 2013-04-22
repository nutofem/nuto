// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/sections/SectionVolume.h"

// create a new section
int NuTo::StructureBase::SectionCreate(const std::string& rType)
{
    // convert section type string to upper case
    std::string SectionTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(SectionTypeString), (int(*)(int)) toupper);

    // get section type from string
    Section::eSectionType SectionType;
    if (SectionTypeString == "TRUSS")
    {
        SectionType = Section::TRUSS;
    }
    else if (SectionTypeString == "PLANE_STRAIN")
    {
        SectionType = Section::PLANE_STRAIN;
    }
    else if (SectionTypeString == "PLANE_STRESS")
    {
        SectionType = Section::PLANE_STRESS;
    }
    else if (SectionTypeString == "VOLUME")
    {
        SectionType = Section::VOLUME;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] invalid section type.");
    }
    return this->SectionCreate(SectionType);
}

// create a new section
int NuTo::StructureBase::SectionCreate(Section::eSectionType rType)
{
	//find unused integer id
    int id(0);
    boost::ptr_map<int,SectionBase>::iterator it = mSectionMap.find(id);
    while (it!=mSectionMap.end())
    {
        id++;
        it = mSectionMap.find(id);
    }

    // create new section
	SectionBase* SectionPtr;
	switch (rType)
	{
	case Section::TRUSS:
		SectionPtr = new SectionTruss();
		break;
	case Section::PLANE_STRAIN:
	case Section::PLANE_STRESS:
		SectionPtr = new SectionPlane(rType);
		break;
	case Section::VOLUME:
		SectionPtr = new SectionVolume();
		break;
	default:
		throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] invalid section type.");
	}

	// add section to map (insert does not allow const keys!!!!)
	this->mSectionMap.insert(id, SectionPtr);
	return id;
}

// delete an existing section
void NuTo::StructureBase::SectionDelete(int rId)
{
    // find section identifier in map
    boost::ptr_map<int,SectionBase>::iterator it = this->mSectionMap.find(rId);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionDelete] Section does not exist.");
    }
    else
    {
        this->mSectionMap.erase(it);
    }
}

// get section pointer from section identifier
NuTo::SectionBase* NuTo::StructureBase::SectionGetSectionPtr(int rId)
{
    boost::ptr_map<int,SectionBase>::iterator it = this->mSectionMap.find(rId);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionGetSectionPtr] Section does not exist.");
    }
    return it->second;
}

// get section pointer from section identifier
const NuTo::SectionBase* NuTo::StructureBase::SectionGetSectionPtr(int rId) const
{
    boost::ptr_map<int,SectionBase>::const_iterator it = this->mSectionMap.find(rId);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionGetSectionPtr] Section does not exist.");
    }
    return it->second;
}

// get section identifier from section pointer
int NuTo::StructureBase::SectionGetId(const NuTo::SectionBase* rSectionPtr) const
{
    for (boost::ptr_map<int,SectionBase>::const_iterator it = mSectionMap.begin(); it!= mSectionMap.end(); it++)
    {
        if (it->second == rSectionPtr)
        {
            return it->first;
        }
    }
    throw MechanicsException("[NuTo::StructureBase::GetSectionId] Section does not exist.");
}

// info routines
void NuTo::StructureBase::SectionInfo(unsigned short rVerboseLevel) const
{
    std::cout << "Number of sections: " << this->GetNumSections() << std::endl;
    for (boost::ptr_map<int,SectionBase>::const_iterator it = mSectionMap.begin(); it!= mSectionMap.end(); it++)
    {
        std::cout << "  Section: " << it->first << std::endl;
        it->second->Info(rVerboseLevel);
    }
}

void NuTo::StructureBase::SectionInfo(int rId, unsigned short rVerboseLevel) const
{
    const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
    std::cout << "  Section: " << rId << std::endl;
    SectionPtr->Info(rVerboseLevel);
}

// set section cross-section area
void NuTo::StructureBase::SectionSetArea(int rId, double rArea)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        SectionPtr->SetArea(rArea);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetArea] error setting section cross-section area.");
        throw e;
    }
}

// get section cross-section area
double NuTo::StructureBase::SectionGetArea(int rId) const
{
    double area;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        area = SectionPtr->GetArea();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionGetArea] error getting section cross-section area.");
        throw e;
    }
    return area;
}

// set section thickness
void NuTo::StructureBase::SectionSetThickness(int rId, double rThickness)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        SectionPtr->SetThickness(rThickness);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetThickness] error setting section thickness.");
        throw e;
    }
}

// get section thickness
double NuTo::StructureBase::SectionGetThickness(int rId) const
{
    double thickness;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        thickness = SectionPtr->GetThickness();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionGetThickness] error getting section thickness.");
        throw e;
    }
    return thickness;
}

// set section dofs
void NuTo::StructureBase::SectionSetDOF(int rId, const std::string& rDOFs)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);

        //now check, if the section is not assigned to any element
        std::vector<ElementBase*> elementVector;
        GetElementsTotal(elementVector);
        for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
        {
        	if (elementVector[countElement]->GetSection()==SectionPtr)
        	{
        		throw MechanicsException("[NuTo::Structure::SectionSetDOF] the section is already assigned to elements, but modification of the properties is only allowed before assigning the section to elements.");
        	}
        }

        // transform string to uppercase
        std::string DOFsUpperCase;
        std::transform(rDOFs.begin(), rDOFs.end(), std::back_inserter(DOFsUpperCase), (int(*)(int)) toupper);

    	// unset all dofs
        SectionPtr->SetIsDisplacementDof(false);
    	SectionPtr->SetIsRotationDof(false);
    	SectionPtr->SetIsTemperatureDof(false);
    	SectionPtr->SetIsNonlocalDamageDof(false);

    	boost::char_separator<char> sep(" ");
        boost::tokenizer< boost::char_separator<char> > tok(DOFsUpperCase, sep);
        for (boost::tokenizer< boost::char_separator<char>  >::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
        {
            if (*beg=="DISPLACEMENTS")
            {
            	SectionPtr->SetIsDisplacementDof(true);
            }
            else if (*beg=="ROTATIONS")
            {
            	SectionPtr->SetIsRotationDof(true);
            }
            else if (*beg=="TEMPERATURE")
            {
            	SectionPtr->SetIsTemperatureDof(true);
            }
            else if (*beg=="NONLOCALDAMAGE")
            {
            	SectionPtr->SetIsNonlocalDamageDof(true);
            }
            else
            {
        		throw MechanicsException("[NuTo::Structure::SectionSetDOF] invalid dof type: " + *beg +".");
            }
        }

    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetDOF] error setting section DOFs.");
        throw e;
    }
}


// set section constitutive inputs
void NuTo::StructureBase::SectionSetInputConstitutive(int rId, const std::string& rDOFs)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);

        //now check, if the section is not assigned to any element
        std::vector<ElementBase*> elementVector;
        GetElementsTotal(elementVector);
        for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
        {
        	if (elementVector[countElement]->GetSection()==SectionPtr)
        	{
        		throw MechanicsException("[NuTo::Structure::SectionSetInputConstitutive] the section is already assigned to elements, but modification of the properties is only allowed before assigning the section to elements.");
        	}
        }

        // transform string to uppercase
        std::string DOFsUpperCase;
        std::transform(rDOFs.begin(), rDOFs.end(), std::back_inserter(DOFsUpperCase), (int(*)(int)) toupper);

    	// unset all dofs
        SectionPtr->SetInputConstitutiveIsTemperature(false);
    	SectionPtr->SetInputConstitutiveIsTemperatureGradient(false);
    	SectionPtr->SetInputConstitutiveIsDeformationGradient(false);

    	boost::char_separator<char> sep(" ");
        boost::tokenizer< boost::char_separator<char> > tok(DOFsUpperCase, sep);
        for (boost::tokenizer< boost::char_separator<char>  >::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
        {
            if (*beg=="TEMPERATURE")
            {
            	SectionPtr->SetInputConstitutiveIsTemperature(true);
            }
            else if (*beg=="TEMPERATUREGRADIENT")
            {
            	SectionPtr->SetInputConstitutiveIsTemperatureGradient(true);
            }
            else if (*beg=="DEFORMATIONGRADIENT")
            {
            	SectionPtr->SetInputConstitutiveIsDeformationGradient(true);
            }
            else if (*beg=="NONLOCALDAMAGE")
            {
            	SectionPtr->SetInputConstitutiveIsDamage(true);
            }
            else
            {
        		throw MechanicsException("[NuTo::Structure::SectionSetInputConstitutive] invalid constitutive input type: " + *beg +".");
            }
        }

    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetInputConstitutive] error setting constitutive input for section.");
        throw e;
    }
}
