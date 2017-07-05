#include "mechanics/interpolationtypes/InterpolationBaseIGA.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

InterpolationBaseIGA::InterpolationBaseIGA(Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, int rDimension)
    : InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
{
}


void InterpolationBaseIGA::Initialize()
{
    mNumNodes = CalculateNumNodes();
    mNumDofs = mNumNodes * Node::GetNumComponents(mDofType, mDimension);
}

