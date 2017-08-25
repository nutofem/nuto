#include <iostream>
#include "BoostUnitTest.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"

class NodeBaseDerived : public NuTo::NodeBase
{
public:
    void SetDofNumber(NuTo::Node::eDof a, int b, int c) override
    {
        std::cout << "Set global DofNumber" << std::endl;
    }

    int GetNumTimeDerivatives(NuTo::Node::eDof) const override
    {
        return 15;
    }

    bool IsDof(NuTo::Node::eDof) const override
    {
        return (true);
    }

    int GetNumDofs() const override
    {
        return 27;
    }

    int GetNum(NuTo::Node::eDof) const override
    {
        return 42;
    }

    int GetDof(NuTo::Node::eDof, int) const override
    {
        return 70;
    }

    const Eigen::VectorXd& Get(NuTo::Node::eDof, int td) const override
    {
        return mValues[td];
    }

    void Set(NuTo::Node::eDof, int td, const Eigen::VectorXd& val) override
    {
        mValues[td] = val;
    }

    std::set<NuTo::Node::eDof> GetDofTypes() const override
    {
        std::set<NuTo::Node::eDof> fakeResult;
        fakeResult.insert(NuTo::Node::eDof::TEMPERATURE);
        return fakeResult;
    }

    void Info(std::ostream& out) const override
    {
    }

    NuTo::NodeBase* Clone() const override
    {
        return new NodeBaseDerived(*this);
    }

    std::vector<Eigen::VectorXd> mValues = {Eigen::Vector2d(1, 1), Eigen::Vector2d(2, 2)};
};

BOOST_AUTO_TEST_CASE(NodeBaseSetGet)
{
    NuTo::NodeBase& nd = *new NodeBaseDerived();

    Eigen::VectorXd expected(1);

    expected(0) = 1.1;
    nd.Set(NuTo::Node::eDof::DISPLACEMENTS, expected(0));
    BOOST_CHECK_EQUAL(nd.Get(NuTo::Node::eDof::DISPLACEMENTS), expected);

    expected(0) = 1.2;
    nd.Set(NuTo::Node::eDof::DISPLACEMENTS, 1, expected(0));
    BOOST_CHECK_EQUAL(nd.Get(NuTo::Node::eDof::DISPLACEMENTS, 1), expected);
}
