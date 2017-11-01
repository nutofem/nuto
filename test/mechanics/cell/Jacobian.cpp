#include "BoostUnitTest.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/elements/ElementShapeFunctions.h"

BOOST_AUTO_TEST_CASE(Jacobian1DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1();
    Eigen::VectorXd coordinates = Eigen::Vector2d({12., 17});
    NuTo::Jacobian jacobian(coordinates, B, 1);
    BOOST_CHECK_CLOSE(jacobian.Det(), 2.5, 1.e-10);
}

BOOST_AUTO_TEST_CASE(Jacobian2DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(Eigen::Vector2d({0, 0}));
    Eigen::VectorXd coordinates = Eigen::VectorXd(8);
    coordinates << 0, 0, 2, 0, 2, 8, 0, 8; // rectangle from (0,0) to (2,8)


    NuTo::Jacobian jacobian(coordinates, B, 2);
    BOOST_CHECK_CLOSE(jacobian.Det(), 4, 1.e-10); // reference area = 4, this area = 16;
}

BOOST_AUTO_TEST_CASE(Jacobian3DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1();
    Eigen::VectorXd coordinates = Eigen::VectorXd(12);
    coordinates << 0, 0, 0, 10, 0, 0, 0, 5, 0, 0, 0, 42;

    NuTo::Jacobian jacobian(coordinates, B, 3);
    BOOST_CHECK_CLOSE(jacobian.Det(), 10 * 5 * 42, 1.e-10);
}

BOOST_AUTO_TEST_CASE(Jacobian1Din2D)
{
    {
        Eigen::MatrixXd B = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1();
        Eigen::VectorXd coordinates = Eigen::VectorXd(4);
        coordinates << 0, 0, 3, 4;

        NuTo::Jacobian jacobian(coordinates, B, 2);
        BOOST_CHECK_CLOSE(jacobian.Det() * 2, 5, 1.e-10);
        //
        // total length of the element is 5 = sqrt(3*3 + 4*4)
        // factor 2 in Det() comes from the integration point weigth
    }
    {
        Eigen::VectorXd coordinates = Eigen::VectorXd(6);
        coordinates << 0, 0, 1.5, 2, 3, 4;

        Eigen::MatrixXd B0 = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(
                Eigen::VectorXd::Constant(1, -std::sqrt(1. / 3.)));
        Eigen::MatrixXd B1 = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(
                Eigen::VectorXd::Constant(1, std::sqrt(1. / 3.)));
        NuTo::Jacobian jacobian0(coordinates, B0, 2);
        NuTo::Jacobian jacobian1(coordinates, B1, 2);
        BOOST_CHECK_CLOSE(jacobian0.Det() + jacobian1.Det(), 5, 1.e-10);
        //
        // total length of the element is 5 = sqrt(3*3 + 4*4)
        // each integration point has weight 1
    }
}

BOOST_AUTO_TEST_CASE(Jacobian2Din3D_TetrahedronLinear)
{
    // define integration points
    std::vector<Eigen::Vector2d> points;
    points.push_back(Eigen::Vector2d(1. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(4. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(1. / 6., 4. / 6.));

    std::vector<Eigen::Vector3d> nodes;
    nodes.push_back(Eigen::Vector3d(1, 0, 0));
    nodes.push_back(Eigen::Vector3d(0, 1, 0));
    nodes.push_back(Eigen::Vector3d(0, 0, 1));
    nodes.push_back(Eigen::Vector3d(1, 1, 1));

    Eigen::VectorXd nodeValuesTetrahedron;
    nodeValuesTetrahedron.resize(12);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
            nodeValuesTetrahedron(3 * i + j) = nodes[i](j);
    }

    for (Eigen::Vector2d& surfaceIP : points)
    {
        for (int surface = 0; surface < 4; surface++)
        {
            int count = 0;
            Eigen::VectorXd nodeValuesTriangle(9);
            Eigen::Vector3d volumetricParameters;
            Eigen::MatrixXd transformation = Eigen::Matrix<double, 3, 2>::Zero();
            switch (surface)
            {
            case 0:
            {
                volumetricParameters << surfaceIP(1), surfaceIP(0), 0; // x,y
                transformation(0, 1) = 1.;
                transformation(1, 0) = 1.;

                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[0](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[1](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[2](i);

                break;
            }
            case 1:
            {
                volumetricParameters << 0, surfaceIP(1), surfaceIP(0); // y, z
                transformation(1, 1) = 1.;
                transformation(2, 0) = 1.;
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[0](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[2](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[3](i);

                break;
            }
            case 2:
            {
                volumetricParameters << surfaceIP(0), 0, surfaceIP(1); // x, z
                transformation(0, 0) = 1.;
                transformation(2, 1) = 1.;
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[0](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[1](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[3](i);

                break;
            }
            case 3:
            {
                volumetricParameters << 1 - surfaceIP(0) - surfaceIP(1), surfaceIP(0), surfaceIP(1); // x, y, z
                transformation(0, 0) = -1.;
                transformation(0, 1) = -1.;
                transformation(1, 0) = 1.;
                transformation(2, 1) = 1.;
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[1](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[2](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[3](i);

                break;
            }
            }

            Eigen::MatrixXd TetrahedronDN = NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1();
            NuTo::Jacobian jacobianTetrahedron(nodeValuesTetrahedron, TetrahedronDN * transformation, 3);

            Eigen::MatrixXd TriangleDN = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(surfaceIP);
            NuTo::Jacobian jacobianTriangle(nodeValuesTriangle, TriangleDN, 3);

            BOOST_CHECK_CLOSE(jacobianTetrahedron.Det(), jacobianTriangle.Det(), 1.e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(Jacobian2Din3D_TetrahedronQuadratic)
{
    // define integration points
    std::vector<Eigen::Vector2d> points;
    points.push_back(Eigen::Vector2d(1. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(4. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(1. / 6., 4. / 6.));

    std::vector<Eigen::Vector3d> nodes;
    nodes.push_back(Eigen::Vector3d(1, 0, 0));
    nodes.push_back(Eigen::Vector3d(0, 1, 0));
    nodes.push_back(Eigen::Vector3d(0, 0, 1));
    nodes.push_back(Eigen::Vector3d(1, 1, 1));
    nodes.push_back(Eigen::Vector3d(0.5, 0.5, 0));
    nodes.push_back(Eigen::Vector3d(0, 0.5, 0.5));
    nodes.push_back(Eigen::Vector3d(0.5, 0, 0.5));
    nodes.push_back(Eigen::Vector3d(1, 0.5, 0.5));
    nodes.push_back(Eigen::Vector3d(0.5, 0.5, 1));
    nodes.push_back(Eigen::Vector3d(0.5, 1, 0.5));

    Eigen::VectorXd nodeValuesTetrahedron;
    nodeValuesTetrahedron.resize(30);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 3; j++)
            nodeValuesTetrahedron(3 * i + j) = nodes[i](j);
    }

    for (Eigen::Vector2d& surfaceIP : points)
    {
        for (int surface = 0; surface < 4; surface++)
        {
            int count = 0;
            Eigen::VectorXd nodeValuesTriangle(18);
            Eigen::Vector3d volumetricParameters;
            Eigen::MatrixXd transformation = Eigen::Matrix<double, 3, 2>::Zero();
            switch (surface)
            {
            case 0:
            {
                volumetricParameters << surfaceIP(1), surfaceIP(0), 0; // x,y
                transformation(0, 1) = 1.;
                transformation(1, 0) = 1.;

                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[0](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[1](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[2](i);

                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[4](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[5](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[6](i);
                break;
            }
            case 1:
            {
                volumetricParameters << 0, surfaceIP(1), surfaceIP(0); // y, z
                transformation(1, 1) = 1.;
                transformation(2, 0) = 1.;
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[0](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[2](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[3](i);

                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[6](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[8](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[7](i);
                break;
            }
            case 2:
            {
                volumetricParameters << surfaceIP(0), 0, surfaceIP(1); // x, z
                transformation(0, 0) = 1.;
                transformation(2, 1) = 1.;
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[0](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[1](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[3](i);

                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[4](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[9](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[7](i);
                break;
            }
            case 3:
            {
                volumetricParameters << 1 - surfaceIP(0) - surfaceIP(1), surfaceIP(0), surfaceIP(1); // x, y, z
                transformation(0, 0) = -1.;
                transformation(0, 1) = -1.;
                transformation(1, 0) = 1.;
                transformation(2, 1) = 1.;
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[1](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[2](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[3](i);

                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[5](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[8](i);
                for (int i = 0; i < 3; i++, count++)
                    nodeValuesTriangle(count) = nodes[9](i);
                break;
            }
            }

            Eigen::MatrixXd TetrahedronDN =
                    NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder2(volumetricParameters);
            NuTo::Jacobian jacobianTetrahedron(nodeValuesTetrahedron, TetrahedronDN * transformation, 3);

            Eigen::MatrixXd TriangleDN = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder2(surfaceIP);
            NuTo::Jacobian jacobianTriangle(nodeValuesTriangle, TriangleDN, 3);


            BOOST_CHECK_CLOSE(jacobianTetrahedron.Det(), jacobianTriangle.Det(), 1.e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(JacobianTransform)
{
    // TODO! How?
}
