#pragma once

using std::cout;
using std::endl;
using NuTo::Constraint::Component;
using NuTo::Constraint::RhsRamp;
using NuTo::Group;
using NuTo::eIntegrationType;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;
using NuTo::eDirection;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
namespace boostFs = boost::filesystem;