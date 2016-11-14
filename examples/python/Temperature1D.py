import numpy as np
import nuto

# Geometry
length = 1.0
# Material
conductivity = 1.0
capacity = 1.0
density = 1.0

def analytic_solution(x, t=0.0):
    return 1.0 + x**2 + 2.0*t

def interpolate(structure, function):
    num_nodes = structure.GetNumNodes()

    coordinates = nuto.DoubleFullVector(1)
    x = np.empty((1, 1))
    for i in range(num_nodes):
        structure.NodeGetCoordinates(i, coordinates)
        coordinates.convrtMatrixToNumpy(x)
        value = function(x[0][0])
        structure.NodeSetTemperature(i, value)
        structure.NodeSetTemperature(i, 1, 2.0)
    return structure

def create_structure(number_of_time_derivatives=0):
    # Geometry/Mesh
    area = 1.0
    number_of_elements = 10

    # create one-dimensional structure
    structure = nuto.Structure(1)
    structure.SetNumTimeDerivatives(number_of_time_derivatives)

    # create section
    section = structure.SectionCreate("Truss")
    structure.SectionSetArea(section, area)

    # create material law
    material = structure.ConstitutiveLawCreate("Heat_Conduction")
    structure.ConstitutiveLawSetParameterDouble(material, "Thermal_Conductivity", conductivity)
    structure.ConstitutiveLawSetParameterDouble(material, "Heat_Capacity", capacity)
    structure.ConstitutiveLawSetParameterDouble(material, "Density", density)

    # create nodes
    node_coordinates = nuto.DoubleFullVector(1)
    for node in range(0, number_of_elements + 1):
        node_coordinates.SetValue(0, 0, node * length/number_of_elements)
        structure.NodeCreate(node, node_coordinates)

    # create interpolation type
    truss_interpolation = structure.InterpolationTypeCreate("Truss1D")
    structure.InterpolationTypeAdd(truss_interpolation, "coordinates", "equidistant1")
    structure.InterpolationTypeAdd(truss_interpolation, "temperature", "equidistant1")
    structure.InterpolationTypeSetIntegrationType(truss_interpolation, "1D2NGauss2Ip")

    # create elements
    element_incidence = nuto.IntFullVector(2)
    for element in range(0, number_of_elements):
        element_incidence.SetValue(0, 0, element)
        element_incidence.SetValue(1, 0, element + 1)
        structure.ElementCreate(truss_interpolation, element_incidence)
        structure.ElementSetSection(element, section)
        structure.ElementSetConstitutiveLaw(element, material)

    structure.ElementTotalConvertToInterpolationType()

    # visualize results
    visualization_group = structure.GroupCreate("Elements")
    structure.GroupAddElementsTotal(visualization_group)

    structure.AddVisualizationComponent(visualization_group, "Temperature")

    return structure


def static_solve(structure):
     # Boundaries
    boundary_temperature = 0.0
    boundary_flux = 10.0

    # set dirichlet boundary condition
    structure.ConstraintLinearSetTemperatureNode(0, boundary_temperature)

    # set Neumann bc
    last_node_id = structure.GetNumNodes() - 1
    direction = nuto.DoubleFullMatrix(1, 1, (1,))
    direction.SetValue(0, 0, 1.0)
    structure.LoadCreateNodeHeatFlux(0, last_node_id, direction, boundary_flux)

    # start analysis
    structure.SolveGlobalSystemStaticElastic(0)
    int_gradient_tmp = structure.BuildGlobalInternalGradient()
    int_gradient = int_gradient_tmp.J.Get("Temperature")
    ext_gradient_tmp = structure.BuildGlobalExternalLoadVector(0)
    ext_gradient = ext_gradient_tmp.J.Get("Temperature")
    residual = nuto.DoubleFullVector(int_gradient - ext_gradient)
    print("Residual: {0}".format(residual.Norm()))

    structure.ExportVtkDataFileElements("Temperature1D_static.vtk")


def transient_solve(structure):

    simulation_time = 1.8
    newmark = nuto.NewmarkDirect(structure)
    newmark.SetTimeStep(.1*simulation_time)
    newmark.SetToleranceForce(1e-4)
    newmark.SetAutomaticTimeStepping(True)

    beta = 2*conductivity/(density*capacity)
    # set dirichlet bc
    boundary_temperature_east = analytic_solution(0.0)
    boundary_temperature_west = analytic_solution(length)
    last_node_id = structure.GetNumNodes() - 1
    bc_west = structure.ConstraintLinearSetTemperatureNode(last_node_id, boundary_temperature_west)
    bc_east = structure.ConstraintLinearSetTemperatureNode(0, boundary_temperature_east)

    end_temp_east = beta*simulation_time + boundary_temperature_east
    temp_east = nuto.DoubleFullMatrix(2, 2)
    temp_east.SetValue(0, 0, 0.0)
    temp_east.SetValue(1, 0, simulation_time)
    temp_east.SetValue(0, 1, boundary_temperature_east)
    temp_east.SetValue(1, 1, end_temp_east)
    newmark.AddTimeDependentConstraint(bc_east, temp_east)

    end_temp_west = beta*simulation_time + boundary_temperature_west
    temp_west = nuto.DoubleFullMatrix(2, 2)
    temp_west.SetValue(0, 0, 0.0)
    temp_west.SetValue(1, 0, simulation_time)
    temp_west.SetValue(0, 1, boundary_temperature_west)
    temp_west.SetValue(1, 1, end_temp_west)
    newmark.AddTimeDependentConstraint(bc_west, temp_west)

    delete_directory = True
    newmark.SetResultDirectory("results_temp_1d", delete_directory)
    newmark.Solve(simulation_time)


def compare_to_analytic(structure):
    num_nodes = structure.GetNumNodes()
    coordinates = nuto.DoubleFullVector(1)
    x = np.empty((1, 1))
    exact_values = np.empty(num_nodes)
    fem_values = np.empty(num_nodes)
    for i in range(num_nodes):
        transient_structure.NodeGetCoordinates(i, coordinates)
        coordinates.convrtMatrixToNumpy(x)
        exact_values[i] = analytic_solution(x[0][0], 1.8)
        fem_values[i] = structure.NodeGetTemperature(i)

    errornorm = np.linalg.norm(exact_values - fem_values, 1) / np.linalg.norm(exact_values, 1)
    print(errornorm)

if __name__ == "__main__":
    static_structure = create_structure()
    static_solve(static_structure)

    transient_structure = create_structure(number_of_time_derivatives=1)
    transient_structure = interpolate(transient_structure, analytic_solution)
    transient_solve(transient_structure)

    compare_to_analytic(transient_structure)
