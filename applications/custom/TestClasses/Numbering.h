#include <vector>
#include <map>

#include <eigen3/Eigen/Core>

#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/nodes/NodeBase.h"

class Numbering
{
public:


    Numbering(){}


    void addGlobalDof(NuTo::Node::eDof rDofType, int rDimension, int rTimeDerivative);


    void addGlobalDof(NuTo::Node::eDof rDofType, int rTimeDerivative, double rDofValue);





    void setGlobalDofValue(NuTo::Node::eDof rDofType, int rGlobalIndex, int rTimeDerivative, double rDofValue);

    void setGlobalDofValueOfNode(NuTo::NodeBase* rNodePtr);

    void setGlobalDofValuesDofType(int rTimeDerivative, NuTo::Node::eDof rDofType, const Eigen::VectorXd& rActiveDofValues,
                            const Eigen::VectorXd& rDependentDofValues);

    int getGlobalIDCount();

    std::vector<int> getGloalIDs(NuTo::Node::eDof rDofType);

    std::vector<double> getGlobalValues(NuTo::Node::eDof rDofType, int rTimeDerivative);



    void printAllDofTypes();

    void printEntries();



    std::vector<std::map<NuTo::Node::eDof, double>> mGlobalIDData;
    //std::vector<std::map<NuTo::Node::eDof, std::vector<double>>>* mGlobalIDData;

    //std::map<NuTo::Node::eDof, std::pair<std::vector<int>, std::map<int, double>>> mGlobalDofData;
    std::map<NuTo::Node::eDof, std::pair<std::vector<int>, std::vector<std::vector<double>>>> mGlobalDofData;


};
