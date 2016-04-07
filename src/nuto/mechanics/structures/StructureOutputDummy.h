/*
 * StructureOutputDummy.h
 *
 *  Created on: 6 Apr 2016
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/structures/StructureOutputBase.h"

namespace NuTo
{

class StructureOutputDummy: public StructureOutputBase
{
public:
    void SetZero() override {}
};

} /* namespace NuTo */
