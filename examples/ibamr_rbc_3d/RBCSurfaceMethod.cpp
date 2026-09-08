#include "RBCSurfaceMethod.h"

#include <ibtk/FEDataManager.h>
#include <ibtk/libmesh_utilities.h>

#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh_base.h>
#include <libmesh/numeric_vector.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace KernelPackRBC
{
using libMesh::DofMap;
using libMesh::EquationSystems;
using libMesh::ExplicitSystem;
using libMesh::MeshBase;
using libMesh::NumericVector;
using libMesh::Point;
using libMesh::VectorValue;

const std::string RBCSurfaceMethod::MEMBRANE_FORCE_SYSTEM_NAME = "MEMBRANE_FORCE_DENSITY";

namespace
{
double
clampUnit(const double x)
{
    return std::max(-1.0, std::min(1.0, x));
}

double
triangleArea(const Point& a, const Point& b, const Point& c)
{
    return 0.5 * ((b - a).cross(c - a)).norm();
}

double
signedTriangleVolume(const Point& a, const Point& b, const Point& c)
{
    return a * b.cross(c) / 6.0;
}
} // namespace

void
RBCSurfaceMethod::initializeMembraneModel(const MembraneParameters& parameters, const unsigned int part)
{
    initializeFEEquationSystems();
    EquationSystems* equation_systems = getFEDataManager(part)->getEquationSystems();
    if (equation_systems->has_system(MEMBRANE_FORCE_SYSTEM_NAME))
    {
        throw std::runtime_error("The membrane force system was initialized more than once.");
    }

    auto& force_system = equation_systems->add_system<ExplicitSystem>(MEMBRANE_FORCE_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        force_system.add_variable("f_membrane_" + std::to_string(d), libMesh::FIRST, libMesh::LAGRANGE);
    }

    std::vector<int> variables(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) variables[d] = static_cast<int>(d);
    std::vector<IBTK::SystemData> system_data(
        1, IBTK::SystemData(MEMBRANE_FORCE_SYSTEM_NAME, variables));
    registerLagSurfaceForceFunction(
        LagSurfaceForceFcnData(interpolateMembraneForce, system_data, nullptr), part);

    d_parameters = parameters;
    d_part = part;
    buildReferenceGeometry(part);
    d_membrane_initialized = true;
}

void
RBCSurfaceMethod::buildReferenceGeometry(const unsigned int part)
{
    EquationSystems* equation_systems = getFEDataManager(part)->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();

    libMesh::dof_id_type max_node_id = 0;
    for (const auto* node : mesh.node_ptr_range()) max_node_id = std::max(max_node_id, node->id());
    d_reference_positions.assign(max_node_id + 1, Point());
    d_reference_lumped_areas.assign(max_node_id + 1, 0.0);
    for (const auto* node : mesh.node_ptr_range()) d_reference_positions[node->id()] = *node;

    d_faces.clear();
    d_reference_area = 0.0;
    d_reference_volume = 0.0;
    std::map<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>,
             std::vector<libMesh::dof_id_type>>
        edge_opposites;

    for (const auto* element : mesh.active_element_ptr_range())
    {
        if (element->n_nodes() != 3)
        {
            throw std::runtime_error("The RBC membrane requires linear triangular elements.");
        }
        Face face{ { element->node_id(0), element->node_id(1), element->node_id(2) } };
        d_faces.push_back(face);

        const Point& a = d_reference_positions[face.node_ids[0]];
        const Point& b = d_reference_positions[face.node_ids[1]];
        const Point& c = d_reference_positions[face.node_ids[2]];
        const double area = triangleArea(a, b, c);
        d_reference_area += area;
        d_reference_volume += signedTriangleVolume(a, b, c);
        for (const auto id : face.node_ids) d_reference_lumped_areas[id] += area / 3.0;

        for (unsigned int e = 0; e < 3; ++e)
        {
            const auto i = face.node_ids[e];
            const auto j = face.node_ids[(e + 1) % 3];
            const auto k = face.node_ids[(e + 2) % 3];
            edge_opposites[std::minmax(i, j)].push_back(k);
        }
    }

    if (d_reference_volume < 0.0) d_reference_volume = -d_reference_volume;
    d_edges.clear();
    d_edges.reserve(edge_opposites.size());
    for (const auto& item : edge_opposites)
    {
        if (item.second.size() != 2)
        {
            throw std::runtime_error("The RBC membrane must be a closed two-manifold surface.");
        }
        const auto i = item.first.first;
        const auto j = item.first.second;
        const auto k = item.second[0];
        const auto l = item.second[1];
        std::array<Point, 4> x{ d_reference_positions[i],
                               d_reference_positions[j],
                               d_reference_positions[k],
                               d_reference_positions[l] };
        d_edges.push_back(
            { i, j, k, l, (x[1] - x[0]).norm(), dihedralAngle(x) });
    }
}

double
RBCSurfaceMethod::dihedralAngle(const std::array<Point, 4>& x)
{
    VectorValue<double> n0 = (x[1] - x[0]).cross(x[2] - x[0]);
    VectorValue<double> n1 = (x[3] - x[0]).cross(x[1] - x[0]);
    const double n0_norm = n0.norm();
    const double n1_norm = n1.norm();
    if (n0_norm <= std::numeric_limits<double>::epsilon() ||
        n1_norm <= std::numeric_limits<double>::epsilon())
    {
        return 0.0;
    }
    n0 /= n0_norm;
    n1 /= n1_norm;
    return std::acos(clampUnit(n0 * n1));
}

void
RBCSurfaceMethod::computeLagrangianForce(const double data_time)
{
    if (!d_membrane_initialized)
    {
        throw std::runtime_error("initializeMembraneModel() must be called before time integration.");
    }
    updateMembraneForce(data_time, d_part);
    IBAMR::IBFESurfaceMethod::computeLagrangianForce(data_time);
}

void
RBCSurfaceMethod::updateMembraneForce(const double /*data_time*/, const unsigned int part)
{
    EquationSystems* equation_systems = getFEDataManager(part)->getEquationSystems();
    auto& coordinate_system = equation_systems->get_system(IBFESurfaceMethod::COORDS_SYSTEM_NAME);
    auto& force_system = equation_systems->get_system<ExplicitSystem>(MEMBRANE_FORCE_SYSTEM_NAME);
    const DofMap& coordinate_dof_map = coordinate_system.get_dof_map();
    const DofMap& force_dof_map = force_system.get_dof_map();

    const NumericVector<double>& coordinate_vector = coordinate_system.get_vector("half");
    std::vector<double> coordinate_values;
    coordinate_vector.localize(coordinate_values);

    std::vector<Point> x(d_reference_positions.size());
    std::vector<unsigned int> dofs;
    for (const auto* node : equation_systems->get_mesh().node_ptr_range())
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            coordinate_dof_map.dof_indices(node, dofs, d);
            x[node->id()](d) = coordinate_values[dofs[0]];
        }
    }

    std::vector<VectorValue<double>> nodal_force(x.size(), VectorValue<double>());

    // In-plane edge elasticity is a discrete, translationally invariant
    // membrane shear model. Each edge retains its reference length.
    for (const Edge& edge : d_edges)
    {
        const VectorValue<double> displacement = x[edge.j] - x[edge.i];
        const double length = displacement.norm();
        if (length <= std::numeric_limits<double>::epsilon()) continue;
        const VectorValue<double> force = d_parameters.shear_modulus *
                                          (length - edge.reference_length) *
                                          displacement / length;
        nodal_force[edge.i] += force;
        nodal_force[edge.j] -= force;
    }

    double current_area = 0.0;
    double signed_volume = 0.0;
    for (const Face& face : d_faces)
    {
        const Point& a = x[face.node_ids[0]];
        const Point& b = x[face.node_ids[1]];
        const Point& c = x[face.node_ids[2]];
        current_area += triangleArea(a, b, c);
        signed_volume += signedTriangleVolume(a, b, c);
    }
    const double orientation = signed_volume < 0.0 ? -1.0 : 1.0;
    const double current_volume = std::abs(signed_volume);
    const double area_factor = d_parameters.area_modulus *
                               (current_area - d_reference_area) /
                               (d_reference_area * d_reference_area);
    const double volume_factor = d_parameters.volume_modulus *
                                 (current_volume - d_reference_volume) /
                                 (d_reference_volume * d_reference_volume);

    // Global area and volume penalties enforce the nearly inextensible,
    // incompressible character of an RBC membrane without tethering its
    // position or orientation in the channel.
    for (const Face& face : d_faces)
    {
        const auto ia = face.node_ids[0];
        const auto ib = face.node_ids[1];
        const auto ic = face.node_ids[2];
        const Point& a = x[ia];
        const Point& b = x[ib];
        const Point& c = x[ic];
        VectorValue<double> normal = (b - a).cross(c - a);
        const double twice_area = normal.norm();
        if (twice_area > std::numeric_limits<double>::epsilon())
        {
            normal /= twice_area;
            nodal_force[ia] -= 0.5 * area_factor * (b - c).cross(normal);
            nodal_force[ib] -= 0.5 * area_factor * (c - a).cross(normal);
            nodal_force[ic] -= 0.5 * area_factor * (a - b).cross(normal);
        }
        nodal_force[ia] -= orientation * volume_factor * b.cross(c) / 6.0;
        nodal_force[ib] -= orientation * volume_factor * c.cross(a) / 6.0;
        nodal_force[ic] -= orientation * volume_factor * a.cross(b) / 6.0;
    }

    // The hinge energy is rotation invariant. A centered local difference is
    // used only for its small four-node gradient, not for the fluid or surface
    // PDE discretizations.
    for (const Edge& edge : d_edges)
    {
        std::array<Point, 4> local{ x[edge.i], x[edge.j], x[edge.k], x[edge.l] };
        const double epsilon = 1.0e-6 * edge.reference_length;
        const std::array<libMesh::dof_id_type, 4> ids{ edge.i, edge.j, edge.k, edge.l };
        for (unsigned int node = 0; node < 4; ++node)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                local[node](d) += epsilon;
                const double angle_plus = dihedralAngle(local);
                const double energy_plus = 0.5 * d_parameters.bending_modulus *
                                           edge.reference_length *
                                           std::pow(angle_plus - edge.reference_dihedral, 2);
                local[node](d) -= 2.0 * epsilon;
                const double angle_minus = dihedralAngle(local);
                const double energy_minus = 0.5 * d_parameters.bending_modulus *
                                            edge.reference_length *
                                            std::pow(angle_minus - edge.reference_dihedral, 2);
                local[node](d) += epsilon;
                nodal_force[ids[node]](d) -= (energy_plus - energy_minus) / (2.0 * epsilon);
            }
        }
    }

    VectorValue<double> residual_force;
    for (const auto* node : equation_systems->get_mesh().node_ptr_range())
        residual_force += nodal_force[node->id()];
    residual_force /= static_cast<double>(equation_systems->get_mesh().n_nodes());

    force_system.solution->zero();
    const auto first_local = force_system.solution->first_local_index();
    const auto last_local = force_system.solution->last_local_index();
    for (const auto* node : equation_systems->get_mesh().node_ptr_range())
    {
        const auto id = node->id();
        const double area = d_reference_lumped_areas[id];
        if (area <= std::numeric_limits<double>::epsilon()) continue;
        const VectorValue<double> force_density = (nodal_force[id] - residual_force) / area;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            force_dof_map.dof_indices(node, dofs, d);
            const auto dof = dofs[0];
            if (first_local <= dof && dof < last_local)
                force_system.solution->set(dof, force_density(d));
        }
    }
    force_system.solution->close();
    force_system.update();
}

void
interpolateMembraneForce(
    VectorValue<double>& force,
    const VectorValue<double>& /*current_normal*/,
    const VectorValue<double>& /*reference_normal*/,
    const libMesh::TensorValue<double>& /*deformation_gradient*/,
    const Point& /*current_position*/,
    const Point& /*reference_position*/,
    libMesh::Elem* /*element*/,
    const unsigned short /*side*/,
    const std::vector<const std::vector<double>*>& variable_data,
    const std::vector<const std::vector<VectorValue<double>>*>& /*variable_gradient_data*/,
    const double /*data_time*/,
    void* /*context*/)
{
    force.zero();
    if (variable_data.empty()) return;
    for (unsigned int d = 0; d < NDIM; ++d) force(d) = (*variable_data[0])[d];
}
} // namespace KernelPackRBC
