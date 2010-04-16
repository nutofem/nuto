// $Id:$

/*******************************************************************************
 Bauhaus-University Weimar
 Author: Joerg F. Unger ,  September 2009
*******************************************************************************/

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "nuto/metamodel/SupportPoints.h"
#include "nuto/metamodel/MetamodelException.h"

using namespace NuTo;

SupportPoints::SupportPoints()
{
    mTransformationBuild = false;
}

SupportPoints::~SupportPoints()
{
    ClearTransformations();
}

//! @brief clear support points
void SupportPoints::Clear()
{
    mSPOrigInput.Resize(0,0);
    mSPTransInput.Resize(0,0);
    mSPTransOutput.Resize(0,0);
    mWeight.resize(0);

    mlTransformationInput.clear();
    mlTransformationOutput.clear();

    mTransformationBuild = false;
}

//! @brief info about support points
void SupportPoints::Info()const
{
    throw MetamodelException("[SupportPoints::Info()] not yet implemented.");
}

void SupportPoints::BuildTransformation()
{
	mSPTransInput = mSPOrigInput;
	for (boost::ptr_list<Transformation>::iterator it = mlTransformationInput.begin(); it!=mlTransformationInput.end();it++)
	{
		it->Build(mSPTransInput);
        it->TransformForward(mSPTransInput);
    }
    
	mSPTransOutput = mSPOrigOutput;
	for (boost::ptr_list<Transformation>::iterator it = mlTransformationOutput.begin(); it!=mlTransformationOutput.end();it++)
	{
		it->Build(mSPTransOutput);
        it->TransformForward(mSPTransOutput);
    }
	mTransformationBuild = true;
}


void SupportPoints::AppendTransformationInput(Transformation* rTransformation)
{
    mlTransformationInput.push_back(rTransformation);
}

void SupportPoints::AppendTransformationOutput(Transformation* rTransformation)
{
    mlTransformationOutput.push_back(rTransformation);
}

void SupportPoints::SetSupportPoints(const FullMatrix<double>& rSPOrigInput, const FullMatrix<double>& rSPOrigOutput)
{
    if (rSPOrigInput.GetNumColumns()!=rSPOrigOutput.GetNumColumns())
    {
        throw MetamodelException("[NuTo::SupportPoints::SetSupportPoints] Number of columns for input and output must be identical (=number of samples).");   
    }
    
    mSPOrigInput = rSPOrigInput;
    mSPOrigOutput  = rSPOrigOutput;

	mWeight.resize(GetNumSupportPoints(),1);
    
    // clear all transformations, since it is no longer for sure that the dimensions remain identical
    mTransformationBuild = false;
    ClearTransformations();
}

//! @brief perform forward transformation for inputs (from orig to transformed)
void SupportPoints::TransformForwardInput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_iterator it=mlTransformationInput.begin();it!=mlTransformationInput.end(); it++)
        it->TransformForward(rCoordinates);
}

//! @brief perform backward transformation for inputs  (from transformed to orig)
void SupportPoints::TransformBackwardInput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_reverse_iterator it=mlTransformationInput.rbegin();it!=mlTransformationInput.rend(); it++)
        it->TransformBackward(rCoordinates);
}

//! @brief perform forward transformation for outputs (from transformed to orig)
//! @brief attention, this is exactly the backwards order, since transformations for outputs are given in revers order
//! @brief orig->trans1->trans2-transformed
void SupportPoints::TransformForwardOutput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_reverse_iterator it=mlTransformationOutput.rbegin();it!=mlTransformationOutput.rend(); it++)
        it->TransformBackward(rCoordinates);
}    

//! @brief perform backward transformation for outputs
//! @brief attention, this is exactly the backwards order, since transformations for given in revers order
void SupportPoints::TransformBackwardOutput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_iterator it=mlTransformationOutput.begin();it!=mlTransformationOutput.end(); it++)
        it->TransformForward(rCoordinates);
}

//! @brief Clears all the transformations for input and output
void SupportPoints::ClearTransformations()
{
    mlTransformationInput.clear();
    mlTransformationOutput.clear();
}
