/*
 * BlockStorageBase.cpp
 *
 *  Created on: 4 Apr 2016
 *      Author: ttitsche
 */

#include "mechanics/dofSubMatrixStorage/BlockStorageBase.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"


NuTo::BlockStorageBase::~BlockStorageBase() {}

int NuTo::BlockStorageBase::GetNumColumns() const
{
    return GetNumColumnsDof(mDofStatus.GetDofTypes());
}

int NuTo::BlockStorageBase::GetNumRows() const
{
    return GetNumRowsDof(mDofStatus.GetDofTypes());
}

int NuTo::BlockStorageBase::GetNumActiveColumns() const
{
    return GetNumColumnsDof(mDofStatus.GetActiveDofTypes());
}

int NuTo::BlockStorageBase::GetNumActiveRows() const
{
    return GetNumRowsDof(mDofStatus.GetActiveDofTypes());
}
