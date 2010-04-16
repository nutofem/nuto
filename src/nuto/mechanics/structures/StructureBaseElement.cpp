#include <assert.h>
#include <boost/tokenizer.hpp>
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief calls ElementCoefficientMatrix_0,
//! renaming only for clarification in mechanical problems for the end user
void NuTo::StructureBase::ElementStiffness(int rElementId, NuTo::FullMatrix<double>& rResult ,
		NuTo::FullMatrix<int>& rGlobalDofsRow,
		NuTo::FullMatrix<int>& rGlobalDofsColumn)const
{
	try
	{
		ElementCoefficientMatrix_0(rElementId, rResult, rGlobalDofsRow, rGlobalDofsColumn);
	}
	catch(NuTo::MechanicsException &e)
	{
		std::stringstream ss;
		ss << rElementId;
		std::string s = ss.str(); //Gets you a C++ STL string
		e.AddMessage("[NuTo::StructureBase::ElementStiffness] Error calculating stiffness of element "
				+ s + ".");
		throw e;
	}
	catch (...)
	{
		std::stringstream ss;
		ss << rElementId;
		std::string s = ss.str(); //Gets you a C++ STL string
		throw MechanicsException("[NuTo::StructureBase::ElementStiffness] Error calculating stiffness of element "
				+ s + ".");
	}
}
//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
void NuTo::StructureBase::ElementCoefficientMatrix_0(int rElementId,
		                 NuTo::FullMatrix<double>& rResult,
		                 NuTo::FullMatrix<int>& rGlobalDofsRow,
		                 NuTo::FullMatrix<int>& rGlobalDofsColumn)const
{
    const ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    std::vector<int> globalDofsRow,
    		         globalDofsColumn;

     try
    {
         bool symmetryFlag;
    	 elementPtr->CalculateCoefficientMatrix_0(rResult, globalDofsRow, globalDofsColumn, symmetryFlag);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
    	e.AddMessage("[NuTo::StructureBase::ElementCoefficientMatrix_0] Error building element matrix for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementCoefficientMatrix_0] Error building element matrix for element " + ss.str() + ".");
    }

    //cast to FullMatrixInt
    rGlobalDofsRow.Resize(globalDofsRow.size(),1);
    memcpy(rGlobalDofsRow.mEigenMatrix.data(),&globalDofsRow[0],globalDofsRow.size()*sizeof(int));

    rGlobalDofsColumn.Resize(globalDofsColumn.size(),1);
    memcpy(rGlobalDofsColumn.mEigenMatrix.data(),&globalDofsRow[0],globalDofsRow.size()*sizeof(int));
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::StructureBase::ElementGradientInternalPotential(int rElementId,
		NuTo::FullMatrix<double>& rResult,
		NuTo::FullMatrix<int>& rGlobalDofsRow)const
{
    const ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    std::vector<int> globalDofsRow;

    try
    {
    	 elementPtr->CalculateGradientInternalPotential(rResult, globalDofsRow);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
    	e.AddMessage("[NuTo::StructureBase::ElementGradientInternalPotential] Error building element vector for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGradientInternalPotential] Error buildung element vector for element " + ss.str() + ".");
    }

    //cast to FullMatrixInt
    rGlobalDofsRow.Resize(globalDofsRow.size(),1);
    memcpy(rGlobalDofsRow.mEigenMatrix.data(),&globalDofsRow[0],globalDofsRow.size()*sizeof(int));

}


//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
void NuTo::StructureBase::ElementCoefficientMatrix_1(int rElementId,
		     NuTo::FullMatrix<double>& rResult,
             NuTo::FullMatrix<int>& rGlobalDofsRow,
             NuTo::FullMatrix<int>& rGlobalDofsColumn)const
{
	throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix_0] To be implemented.");
}

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
void NuTo::StructureBase::ElementCoefficientMatrix_2(int rElementId,
		     NuTo::FullMatrix<double>& rResult,
             NuTo::FullMatrix<int>& rGlobalDofsRow,
             NuTo::FullMatrix<int>& rGlobalDofsColumn)const
{
	throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix_0] To be implemented.");
}

//! @brief sets the constitutive law of a single element
//! @param rElementIdent identifier for the element
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementSetConstitutiveLaw(int rElementId, const std::string& rConstitutiveLawIdent)
{
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<std::string,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    try
    {
    	ElementSetConstitutiveLaw(elementPtr,itConstitutive->second);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element " + ss.str() + ".");
    }
}

//! @brief sets the constitutive law of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementGroupSetConstitutiveLaw(const std::string& rGroupIdent, const std::string& rConstitutiveLawIdent)
{
	boost::ptr_map<std::string,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::GroupBase::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<std::string,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetConstitutiveLaw(*itElement,itConstitutive->second);
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(*itElement);
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting constitutive law  for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(*itElement);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting constitutive law for element " + ss.str() + ".");
        }
    }
}

//! @brief sets the constitutive law of a all elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementTotalSetConstitutiveLaw(const std::string& rConstitutiveLawIdent)
{
    boost::ptr_map<std::string,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetConstitutiveLaw(elementVector[countElement],itConstitutive->second);
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting constitutive law  for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting constitutive law for element "
        			   + ss.str() + ".");
        }
    }
}

//! @brief sets the constitutive law of a single element
//! @param rElement element pointer
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive)
{
    rElement->SetConstitutiveLaw(rConstitutive);
}


//! @brief sets the section of a single element
//! @param rElementIdent identifier for the element
//! @param rConstitutiveLawIdent identifier for the section
void NuTo::StructureBase::ElementSetSection(int rElementId, const std::string& rSectionIdent)
{
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<std::string,SectionBase>::iterator itSection = mSectionMap.find(rSectionIdent);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetSection] Section with the given identifier does not exist.");

    try
    {
    	ElementSetSection(elementPtr,itSection->second);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetSection] Error setting section for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetSection] Error setting section for element " + ss.str() + ".");
    }
}

//! @brief sets the section of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementGroupSetSection(const std::string& rGroupIdent, const std::string& rSectionIdent)
{
	boost::ptr_map<std::string,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::GroupBase::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<std::string,SectionBase>::iterator itSection = mSectionMap.find(rSectionIdent);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Section with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetSection(*itElement,itSection->second);
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(*itElement);
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting section for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(*itElement);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting section for element " + ss.str() + ".");
        }
    }
}

//! @brief sets the section for all elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementTotalSetSection(const std::string& rSectionIdent)
{
    boost::ptr_map<std::string,SectionBase>::iterator itSection = mSectionMap.find(rSectionIdent);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Section with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetSection(elementVector[countElement],itSection->second);
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting section  for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting section for element "
        			   + ss.str() + ".");
        }
    }
}

//! @brief modifies the material of a single element
//! @param rElement element pointer
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetSection(ElementBase* rElement, SectionBase* rSection)
{
    rElement->SetSection(rSection);
}


//! @brief modifies the section of a single element
//! @param rElementIdent identifier for the element
//! @param rSectionIdent identifier for the section
void NuTo::StructureBase::ElementSetIntegrationType(int rElementId,
		const std::string& rIntegrationTypeIdent)
{
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    try
    {
    	ElementSetIntegrationType(elementPtr,GetPtrIntegrationType(rIntegrationTypeIdent));
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetIntegrationType] Error setting integration type for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetIntegrationType] Error setting integration type for element " + ss.str() + ".");
    }

}

//! @brief modifies the integration type of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rIntegrationTypeIdent identifier for the integration type
void NuTo::StructureBase::ElementGroupSetIntegrationType(const std::string& rGroupIdent,
		const std::string& rIntegrationTypeIdent)
{
	boost::ptr_map<std::string,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetIntegrationType] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::GroupBase::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetIntegrationType] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetIntegrationType(*itElement,GetPtrIntegrationType(rIntegrationTypeIdent));
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(*itElement);
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetIntegrationType] Error setting integration type for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(*itElement);
       	    throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetIntegrationType] Error setting integration type for element " + ss.str() + ".");
        }
    }

}

//! @brief modifies the section of a all elements
//! @param rSectionIdent identifier for the section
void NuTo::StructureBase::ElementTotalSetIntegrationType(const std::string& rIntegrationTypeIdent)
{
    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetIntegrationType(elementVector[countElement],GetPtrIntegrationType(rIntegrationTypeIdent));
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetIntegrationType] Error setting integration type for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetIntegrationType] Error setting integration type for element "
        			   + ss.str() + ".");
        }
    }
}

//! @brief modifies the integration type of a single element
//! @param rElement element pointer
//! @param rIntegrationType integration type
void NuTo::StructureBase::ElementSetIntegrationType(ElementBase* rElement, const IntegrationTypeBase* rIntegrationType)
{
	rElement->SetIntegrationType(rIntegrationType);
}

//! @brief calculates the engineering strain
//! @param rElemIdent  identifier for the element
//! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
void NuTo::StructureBase::ElementGetEngineeringStrain(int rElementId, FullMatrix<double>& rEngineeringStrain)const
{
    const ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    try
    {
    	elementPtr->GetEngineeringStrain(rEngineeringStrain);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetEngineeringStrain] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] Error getting engineering strain for element " + ss.str() + ".");
    }
}

//! @brief calculates the engineering plastic strain
//! @param rElemIdent  identifier for the element
//! @param rEngineerungStrain engineering plastic strain (return value, always 6xnumIp matrix)
void NuTo::StructureBase::ElementGetEngineeringPlasticStrain(int rElementId, FullMatrix<double>& rEngineeringStrain)const
{
    const ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    try
    {
    	elementPtr->GetEngineeringPlasticStrain(rEngineeringStrain);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] Error getting engineering strain for element " + ss.str() + ".");
    }
}

//! @brief calculates the engineering stress
//! @param rElemIdent  identifier for the element
//! @param rEngineeringStress engineering stress (return value, always 6xnumIp matrix)
void NuTo::StructureBase::ElementGetEngineeringStress(int rElementId, FullMatrix<double>& rEngineeringStress)const
{
    const ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
    	elementPtr->GetEngineeringStress(rEngineeringStress);
    }
    catch(NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetEngineeringStrain] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetIntegrationType] Error getting engineering strain for element " + ss.str() + ".");
    }
}

//! @brief updates the history data of a all elements
void NuTo::StructureBase::ElementTotalUpdateStaticData()
{
    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	elementVector[countElement]->UpdateStaticData();
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalUpdateStaticData] Error updating static data for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalUpdateStaticData] Error updating static data for element "
        			   + ss.str() + ".");
        }
    }
}
