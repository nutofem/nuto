#include "Numbering.h"
#include <iostream>
#include <vector>
#include <map>

#include <eigen3/Eigen/Core>

#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/nodes/NodeBase.h"



void Numbering::addGlobalDof(NuTo::Node::eDof rDofType, int rDimension, int rTimeDerivative)
{
    rTimeDerivative++;
    auto dofIDsAndValues = mGlobalDofData.find(rDofType);
    if (dofIDsAndValues == mGlobalDofData.end())
    {
        std::vector<double> newValuesPerTimeDerivative(rDimension);
        std::vector<std::vector<double>> newValues(rTimeDerivative);
        for (int i = 0; i < rTimeDerivative; ++i)
        {
            newValues[i] = newValuesPerTimeDerivative;
        }

        std::vector<int> newIndices(rDimension);

        std::pair<std::vector<int>, std::vector<std::vector<double>>> newEntry;
        newEntry.first = newIndices;
        newEntry.second = newValues;

        mGlobalDofData.insert(std::pair<NuTo::Node::eDof, std::pair<std::vector<int>, std::vector<std::vector<double>>>>(rDofType, newEntry));
    }
    else
    {
        std::vector<int> actualIDs = dofIDsAndValues->second.first;
        std::vector<std::vector<double>> actualValues = dofIDsAndValues->second.second;
        std::vector<int>::iterator maxIt;
        maxIt = std::max_element(actualIDs.begin() ,actualIDs.end());
        int actualMaxID = *maxIt;

        actualIDs.push_back(actualMaxID + 1);

        for (int i = 0; i < rTimeDerivative; ++i)
        {
            actualValues[i].push_back(0);
        }
    }
}

void Numbering::addGlobalDof(NuTo::Node::eDof rDofType, int rTimeDerivative, double rDofValue)
{
    //mGlobalIDData.push_back(rDofValuesMap);

}

void Numbering::setGlobalDofValue(NuTo::Node::eDof rDofType, int rGlobalIndex, int rTimeDerivative, double rDofValue)
{
    auto dofIDsAndValues = mGlobalDofData.find(rDofType);
    if (dofIDsAndValues == mGlobalDofData.end())
    {
        //No entry with given doftype found
    }
    else
    {
        std::vector<int> actualIDs = dofIDsAndValues->second.first;
        std::vector<std::vector<double>> actualValues = dofIDsAndValues->second.second;

        auto foundIndex = std::find(actualIDs.begin(), actualIDs.end(), rGlobalIndex);
        if (foundIndex == actualIDs.end())
        {
            //No entry with given index found
        }
        else
        {
            int position = std::distance(actualIDs.begin(), foundIndex);
            actualValues[rTimeDerivative][position] = rDofValue;
        }
    }
}



void Numbering::setGlobalDofValueOfNode(NuTo::NodeBase *rNodePtr)
//void Numbering::setGlobalDofValueOfNode(NuTo::NodeDof* rNodePtr)
{
    Eigen::VectorXd dofValues;

    for (NuTo::Node::eDof dofType : rNodePtr->GetDofTypes())
    {
        dofValues = rNodePtr->Get(dofType);


        for (int i = 0; i < dofValues.rows(); ++i)
        {

        }
    }
}

void Numbering::setGlobalDofValuesDofType(int rTimeDerivative, NuTo::Node::eDof rDofType, const Eigen::VectorXd &rActiveDofValues, const Eigen::VectorXd &rDependentDofValues)
{

}

int Numbering::getGlobalIDCount()
{
    return int(mGlobalDofData.size());
}

std::vector<int> Numbering::getGloalIDs(NuTo::Node::eDof rDofType)
{
    auto dofIDsAndValues = mGlobalDofData.find(rDofType);
    if (dofIDsAndValues == mGlobalDofData.end())
    {
        //No entry with given dof type
    }
    else
    {
        return dofIDsAndValues->second.first;
    }
}

std::vector<double> Numbering::getGlobalValues(NuTo::Node::eDof rDofType, int rTimeDerivative)
{
    auto dofIDsAndValues = mGlobalDofData.find(rDofType);
    if (dofIDsAndValues == mGlobalDofData.end())
    {
        //No entry with given dof type
    }
    else
    {
        std::vector<std::vector<double>> actualValues = dofIDsAndValues->second.second;
        return actualValues[rTimeDerivative];
    }
}

void Numbering::printAllDofTypes()
{
    for(auto it : mGlobalDofData) {
        std::cout << NuTo::Node::DofToString(it.first) << "\n";
    }
}

void Numbering::printEntries()
{

}

