#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"

namespace NuTo
{

class Element
{
public:
    Element(std::vector<NuTo::NodeSimple*> rNodes)
        : mNodes(rNodes)
    {
    }

    Eigen::VectorXd ExtractNodeValues() const
    {
        int dim = mNodes[0]->GetNumValues();
        Eigen::VectorXd nodeValues(mNodes.size() * dim);
        for (size_t i = 0; i < mNodes.size(); ++i)
            nodeValues.segment(dim * i, dim) = mNodes[i]->GetValues();
        return nodeValues;
    }

private:
    std::vector<NuTo::NodeSimple*> mNodes;
};
} /* NuTo */ 
