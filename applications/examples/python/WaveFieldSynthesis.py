import nuto
import numpy as np
from functools import partial


def smooth_dirac_constraint_function(time, peak_time):
    a = 0.025
    return 1. * np.exp(- ((time - peak_time) / a) ** 2)


def line_distance(l_0, n, p):
    """Calculate the distance from line l_0 + t * n to the point p."""
    v = l_0 - p
    return np.linalg.norm(v - v.dot(n) * n)


LENGTH_X = 3.
LENGTH_Y = 3.
LENGTH_Z = 1.
ELEMENT_SIZE = 0.20
# ELEMENT_SIZE = 0.025


def setup_geometry(s):
    s.SetNumTimeDerivatives(2)
    s.SetShowTime(False)
    s.SetVerboseLevel(10)

    mesh_info = nuto.MeshGenerator_Grid(s, [LENGTH_X, LENGTH_Y], [int(LENGTH_X / ELEMENT_SIZE), int(LENGTH_Y / ELEMENT_SIZE) ])
    # s.InterpolationTypeAdd(mesh_info[1], "ELECTRICPOTENTIAL", "EQUIDISTANT2")
    s.InterpolationTypeAdd(mesh_info[1], "ELECTRICPOTENTIAL", "LOBATTO4")
    s.ElementTotalConvertToInterpolationType()

    law_id = s.ConstitutiveLawCreate("LINEAR_DIELECTRIC")
    dielectric_tensor = np.eye(3)
    s.ConstitutiveLawSetParameterMatrixDouble(law_id, "DIELECTRIC_TENSOR", -dielectric_tensor)
    s.ElementTotalSetConstitutiveLaw(law_id)

    s.ElementTotalSetSection(nuto.SectionPlane_Create(LENGTH_Z, False))

    all_elements = s.GroupGetElementsTotal()
    s.AddVisualizationComponent(all_elements, "ElectricPotential")


def run_time_integration(s, result_directory, end_time):
    rk4 = nuto.RungeKutta4(s)
    ti = nuto.RungeKutta4(s)

    time_step = rk4.CalculateCriticalTimeStep()
    print("Critical time step", time_step)
    ti.SetTimeStep(time_step / 2.)

    ti.PostProcessing().SetResultDirectory(result_directory, True)
    ti.Solve(end_time)


def define_bc_point_source():
    """Define the boundary conditions for a line source.

    This currently assumes that the point source only influences the WEST
    border of the domain. EAST, NORTH and SOUTH are in the _shadow_ of WEST
    """
    nodes_west = s.GroupGetNodesAtCoordinate(nuto.eDirection_X, 0)

    source_position = np.array([-1.5, 1.5])

    for node_id in nodes_west.GetMemberIds():
        node = s.NodeGetNodePtr(node_id)
        pos = node.Get(nuto.eDof_COORDINATES).flatten()
        distance_to_source = np.linalg.norm(pos - source_position)
        time_to_source = distance_to_source - 1.5  # assumens: wave speed = 1
        s.Constraints().Add(nuto.eDof_ELECTRICPOTENTIAL, nuto.Value(node, partial(smooth_dirac_constraint_function, time_to_source)))


def define_bc_line_source(l_0, n):
    """Define the boundary conditions for a line source.

    Wave coming from the line defined by point l_0 and normal vector n.
    Assumes that l_0 is in the SOUTH_WEST of the domain and n pointing towards
    it. So the bcs are applied to the SOUTH and the WEST boundary and NORTH and
    EAST are in their shadow
    """
    nodes_west = s.GroupGetNodesAtCoordinate(nuto.eDirection_X, 0)
    nodes_south = s.GroupGetNodesAtCoordinate(nuto.eDirection_Y, 0)
    id_west = s.GroupGetId(nodes_west)
    id_south = s.GroupGetId(nodes_south)
    id_union = s.GroupUnion(id_west, id_south)

    for node_id in s.GroupGetMemberIds(id_union):
        node = s.NodeGetNodePtr(node_id)
        pos = node.Get(nuto.eDof_COORDINATES).flatten()

        distance_to_source = line_distance(l_0, n, pos)
        time_to_source = distance_to_source  # assumes: wave speed = 1
        s.Constraints().Add(nuto.eDof_ELECTRICPOTENTIAL, nuto.Value(node, partial(smooth_dirac_constraint_function, time_to_source)))


if __name__ == "__main__":
    s = nuto.Structure(2)
    setup_geometry(s)

    define_bc_point_source()

    run_time_integration(s, "./WaveFieldSynthesisResults", .5)
