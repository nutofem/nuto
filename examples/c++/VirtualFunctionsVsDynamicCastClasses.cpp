// $Id$ 
// VirtualFunctionsVsDynamicCastClasses.cpp
// created May 10, 2010 by Joerg F. Unger

#include "VirtualFunctionsVsDynamicCastClasses.h"

void BaseA::FunctionA(NuToVec2 rValue)
{
    rValue[0] = mA[0];
    rValue[1] = mA[1];
}

void BaseA::FunctionB(NuToVec2 rValue)
{
	throw 1;
}

BaseA* BaseA::asA()
{
	return this;
}

void BaseB::FunctionB(NuToVec3 rValue)
{
    rValue[0] = mB[0];
    rValue[1] = mB[1];
    rValue[1] = mB[2];
}

BaseB* BaseB::asB()
{
	return this;
}

AB* AB::asAB()
{
	return this;
}

void AB::FunctionB(NuToVec3 rValue)
{
    rValue[0] = mB[0];
    rValue[1] = mB[1];
    rValue[2] = mB[2];
}
