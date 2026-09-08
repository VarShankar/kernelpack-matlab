#include "RBCSurfaceMethod.h"

#include <SAMRAI_config.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/replicated_mesh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
using Face = std::array<unsigned int, 3>;

struct SurfaceMeshData
{
    std::vector<libMesh::Point> vertices;
    std::vector<Face> faces;
};

unsigned int
midpointVertex(const unsigned int i,
               const unsigned int j,
               std::vector<libMesh::Point>& vertices,
               std::map<std::pair<unsigned int, unsigned int>, unsigned int>& cache)
{
    const auto key = std::minmax(i, j);
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    libMesh::Point midpoint = 0.5 * (vertices[i] + vertices[j]);
    midpoint /= midpoint.norm();
    const unsigned int id = static_cast<unsigned int>(vertices.size());
    vertices.push_back(midpoint);
    cache[key] = id;
    return id;
}

SurfaceMeshData
makeIcosphere(const unsigned int subdivisions)
{
    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    SurfaceMeshData mesh{
        { { -1, phi, 0 }, { 1, phi, 0 },  { -1, -phi, 0 }, { 1, -phi, 0 },
          { 0, -1, phi }, { 0, 1, phi },  { 0, -1, -phi }, { 0, 1, -phi },
          { phi, 0, -1 }, { phi, 0, 1 },  { -phi, 0, -1 }, { -phi, 0, 1 } },
        { { 0, 11, 5 }, { 0, 5, 1 },  { 0, 1, 7 },   { 0, 7, 10 }, { 0, 10, 11 },
          { 1, 5, 9 },  { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
          { 3, 9, 4 },  { 3, 4, 2 },  { 3, 2, 6 },   { 3, 6, 8 },  { 3, 8, 9 },
          { 4, 9, 5 },  { 2, 4, 11 }, { 6, 2, 10 },  { 8, 6, 7 },  { 9, 8, 1 } } };

    for (auto& vertex : mesh.vertices) vertex /= vertex.norm();
    for (unsigned int level = 0; level < subdivisions; ++level)
    {
        std::map<std::pair<unsigned int, unsigned int>, unsigned int> cache;
        std::vector<Face> refined_faces;
        refined_faces.reserve(4 * mesh.faces.size());
        for (const Face& face : mesh.faces)
        {
            const unsigned int a = midpointVertex(face[0], face[1], mesh.vertices, cache);
            const unsigned int b = midpointVertex(face[1], face[2], mesh.vertices, cache);
            const unsigned int c = midpointVertex(face[2], face[0], mesh.vertices, cache);
            refined_faces.push_back({ face[0], a, c });
            refined_faces.push_back({ face[1], b, a });
            refined_faces.push_back({ face[2], c, b });
            refined_faces.push_back({ a, b, c });
        }
        mesh.faces.swap(refined_faces);
    }
    return mesh;
}

libMesh::Point
mapToBiconcaveRBC(const libMesh::Point& q,
                  const libMesh::Point& center,
                  const double radius)
{
    // Evans-Fung biconcave profile. The symmetry axis is streamwise so the
    // membrane begins as a disc whose broad face is normal to the flow.
    constexpr double c0 = 0.207;
    constexpr double c1 = 2.002;
    constexpr double c2 = -1.123;
    const double radial_fraction = std::sqrt(q(0) * q(0) + q(1) * q(1));
    const double radial_squared = radial_fraction * radial_fraction;
    const double half_thickness = 0.5 * radius * std::sqrt(std::max(0.0, 1.0 - radial_squared)) *
                                  (c0 + c1 * radial_squared + c2 * radial_squared * radial_squared);
    const double axial = q(2) < 0.0 ? -half_thickness : half_thickness;
    return center + libMesh::Point(axial, radius * q(0), radius * q(1));
}

std::vector<libMesh::Point>
buildRBCMesh(libMesh::ReplicatedMesh& mesh,
             const unsigned int subdivisions,
             const libMesh::Point& center,
             const double radius)
{
    SurfaceMeshData surface = makeIcosphere(subdivisions);
    for (auto& vertex : surface.vertices) vertex = mapToBiconcaveRBC(vertex, center, radius);

    double signed_volume = 0.0;
    for (const Face& face : surface.faces)
    {
        signed_volume += surface.vertices[face[0]] *
                         surface.vertices[face[1]].cross(surface.vertices[face[2]]) / 6.0;
    }
    if (signed_volume < 0.0)
        for (Face& face : surface.faces) std::swap(face[1], face[2]);

    for (unsigned int id = 0; id < surface.vertices.size(); ++id)
        mesh.add_point(surface.vertices[id], id);
    for (const Face& face : surface.faces)
    {
        std::unique_ptr<libMesh::Elem> element = libMesh::Elem::build(libMesh::TRI3);
        for (unsigned int i = 0; i < 3; ++i) element->set_node(i) = mesh.node_ptr(face[i]);
        mesh.add_elem(element.release());
    }
    mesh.prepare_for_use();
    std::vector<libMesh::Point> material_coordinates(mesh.max_node_id());
    for (const auto* node : mesh.node_ptr_range())
    {
        const double q0 = ((*node)(1) - center(1)) / radius;
        const double q1 = ((*node)(2) - center(2)) / radius;
        const double radial_squared = std::min(1.0, q0 * q0 + q1 * q1);
        const double q2_magnitude = std::sqrt(std::max(0.0, 1.0 - radial_squared));
        const double axial_displacement = (*node)(0) - center(0);
        const double q2 = axial_displacement < 0.0 ? -q2_magnitude : q2_magnitude;
        libMesh::Point material(q0, q1, q2);
        material /= material.norm();
        const double cy = std::cos(0.37);
        const double sy = std::sin(0.37);
        const double cz = std::cos(0.23);
        const double sz = std::sin(0.23);
        const libMesh::Point rotated_y(cy * material(0) + sy * material(2),
                                       material(1),
                                       -sy * material(0) + cy * material(2));
        material_coordinates[node->id()] = libMesh::Point(
            cz * rotated_y(0) - sz * rotated_y(1),
            sz * rotated_y(0) + cz * rotated_y(1),
            rotated_y(2));
    }
    return material_coordinates;
}

std::vector<libMesh::Point>
localizeVectorSystem(libMesh::EquationSystems& equation_systems, const std::string& system_name)
{
    auto& system = equation_systems.get_system(system_name);
    const libMesh::DofMap& dof_map = system.get_dof_map();
    std::vector<double> values;
    system.solution->localize(values);
    std::vector<libMesh::Point> result(equation_systems.get_mesh().max_node_id());
    std::vector<unsigned int> dofs;
    for (const auto* node : equation_systems.get_mesh().node_ptr_range())
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(node, dofs, d);
            result[node->id()](d) = values[dofs[0]];
        }
    }
    return result;
}

void
writeConnectivity(const libMesh::MeshBase& mesh,
                  const std::vector<libMesh::Point>& material_coordinates,
                  const std::filesystem::path& output_directory)
{
    if (IBTK::IBTK_MPI::getRank() != 0) return;
    std::ofstream faces(output_directory / "faces.csv");
    faces << "node_1,node_2,node_3\n";
    for (const auto* element : mesh.active_element_ptr_range())
        faces << element->node_id(0) << ',' << element->node_id(1) << ',' << element->node_id(2) << '\n';

    std::ofstream reference(output_directory / "reference.csv");
    reference << std::setprecision(17) << "node_id,x,y,z\n";
    for (const auto* node : mesh.node_ptr_range())
        reference << node->id() << ',' << (*node)(0) << ',' << (*node)(1) << ',' << (*node)(2) << '\n';

    std::ofstream material(output_directory / "material.csv");
    material << std::setprecision(17) << "node_id,u_x,u_y,u_z\n";
    for (const auto* node : mesh.node_ptr_range())
    {
        const libMesh::Point& u = material_coordinates[node->id()];
        material << node->id() << ',' << u(0) << ',' << u(1) << ',' << u(2) << '\n';
    }

    std::ofstream diagnostics(output_directory / "diagnostics.csv");
    diagnostics << "step,time,area,relative_area_change,volume,relative_volume_change,minimum_edge,maximum_edge,centroid_x,centroid_y,centroid_z\n";
}

void
writeTrajectoryFrame(libMesh::EquationSystems& equation_systems,
                     const int step,
                     const double time,
                     const std::filesystem::path& output_directory,
                     double& reference_area,
                     double& reference_volume)
{
    const std::vector<libMesh::Point> positions = localizeVectorSystem(
        equation_systems, IBAMR::IBFESurfaceMethod::COORDS_SYSTEM_NAME);
    const std::vector<libMesh::Point> velocities = localizeVectorSystem(
        equation_systems, IBAMR::IBFESurfaceMethod::VELOCITY_SYSTEM_NAME);
    if (IBTK::IBTK_MPI::getRank() != 0) return;

    std::ostringstream filename;
    filename << "frame_" << std::setw(6) << std::setfill('0') << step << ".csv";
    std::ofstream frame(output_directory / filename.str());
    frame << std::setprecision(17) << "node_id,x,y,z,u,v,w\n";
    libMesh::Point centroid;
    for (const auto* node : equation_systems.get_mesh().node_ptr_range())
    {
        const auto id = node->id();
        centroid += positions[id];
        frame << id << ',' << positions[id](0) << ',' << positions[id](1) << ',' << positions[id](2) << ','
              << velocities[id](0) << ',' << velocities[id](1) << ',' << velocities[id](2) << '\n';
    }
    centroid /= static_cast<double>(equation_systems.get_mesh().n_nodes());

    double area = 0.0;
    double signed_volume = 0.0;
    double minimum_edge = std::numeric_limits<double>::max();
    double maximum_edge = 0.0;
    std::set<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>> edges;
    for (const auto* element : equation_systems.get_mesh().active_element_ptr_range())
    {
        const auto i = element->node_id(0);
        const auto j = element->node_id(1);
        const auto k = element->node_id(2);
        area += 0.5 * ((positions[j] - positions[i]).cross(positions[k] - positions[i])).norm();
        signed_volume += positions[i] * positions[j].cross(positions[k]) / 6.0;
        edges.insert(std::minmax(i, j));
        edges.insert(std::minmax(j, k));
        edges.insert(std::minmax(k, i));
    }
    for (const auto& edge : edges)
    {
        const double length = (positions[edge.second] - positions[edge.first]).norm();
        minimum_edge = std::min(minimum_edge, length);
        maximum_edge = std::max(maximum_edge, length);
    }
    const double volume = std::abs(signed_volume);
    if (reference_area == 0.0)
    {
        reference_area = area;
        reference_volume = volume;
    }
    std::ofstream diagnostics(output_directory / "diagnostics.csv", std::ios::app);
    diagnostics << std::setprecision(17) << step << ',' << time << ',' << area << ','
                << (area - reference_area) / reference_area << ',' << volume << ','
                << (volume - reference_volume) / reference_volume << ',' << minimum_edge << ',' << maximum_edge
                << ',' << centroid(0) << ',' << centroid(1) << ',' << centroid(2) << '\n';
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const libMesh::LibMeshInit& libmesh_init = ibtk_init.getLibMeshInit();

    {
        Pointer<IBTK::AppInitializer> app_initializer = new IBTK::AppInitializer(argc, argv, "rbc3d.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        libMesh::ReplicatedMesh membrane_mesh(libmesh_init.comm(), 2);
        const unsigned int subdivisions = input_db->getIntegerWithDefault("RBC_SUBDIVISIONS", 3);
        const libMesh::Point center(input_db->getDoubleWithDefault("RBC_CENTER_X", 0.75),
                                    input_db->getDoubleWithDefault("RBC_CENTER_Y", 0.35),
                                    input_db->getDoubleWithDefault("RBC_CENTER_Z", 0.50));
        const std::vector<libMesh::Point> material_coordinates = buildRBCMesh(
            membrane_mesh,
            subdivisions,
            center,
            input_db->getDoubleWithDefault("RBC_RADIUS", 0.18));

        Pointer<IBAMR::INSStaggeredHierarchyIntegrator> fluid_integrator =
            new IBAMR::INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<KernelPackRBC::RBCSurfaceMethod> membrane_method =
            new KernelPackRBC::RBCSurfaceMethod(
                "IBFESurfaceMethod",
                app_initializer->getComponentDatabase("IBFESurfaceMethod"),
                &membrane_mesh,
                app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                true,
                app_initializer->getRestartReadDirectory(),
                app_initializer->getRestartRestoreNumber());

        KernelPackRBC::MembraneParameters membrane_parameters;
        membrane_parameters.shear_modulus = input_db->getDouble("RBC_SHEAR_MODULUS");
        membrane_parameters.bending_modulus = input_db->getDouble("RBC_BENDING_MODULUS");
        membrane_parameters.area_modulus = input_db->getDouble("RBC_AREA_MODULUS");
        membrane_parameters.volume_modulus = input_db->getDouble("RBC_VOLUME_MODULUS");
        membrane_method->initializeMembraneModel(membrane_parameters);

        Pointer<IBAMR::IBHierarchyIntegrator> time_integrator =
            new IBAMR::IBExplicitHierarchyIntegrator(
                "IBHierarchyIntegrator",
                app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                membrane_method,
                fluid_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize",
            time_integrator,
            app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm",
            app_initializer->getComponentDatabase("GriddingAlgorithm"),
            error_detector,
            box_generator,
            load_balancer);

        Pointer<IBTK::CartGridFunction> initial_velocity = new IBTK::muParserCartGridFunction(
            "initial_velocity",
            app_initializer->getComponentDatabase("VelocityInitialConditions"),
            grid_geometry);
        fluid_integrator->registerVelocityInitialConditions(initial_velocity);

        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> velocity_boundary_conditions(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            velocity_boundary_conditions[d] = new IBTK::muParserRobinBcCoefs(
                "velocity_bc_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry);
        }
        fluid_integrator->registerPhysicalBoundaryConditions(velocity_boundary_conditions);

        Pointer<IBTK::CartGridFunction> body_force = new IBTK::muParserCartGridFunction(
            "body_force", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
        time_integrator->registerBodyForceFunction(body_force);

        const bool dump_visualization = app_initializer->dumpVizData();
        const int visualization_interval = app_initializer->getVizDumpInterval();
        Pointer<VisItDataWriter<NDIM>> visit_writer = app_initializer->getVisItDataWriter();
        if (dump_visualization && visit_writer) time_integrator->registerVisItDataWriter(visit_writer);

        membrane_method->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        libMesh::EquationSystems* equation_systems = membrane_method->getFEDataManager()->getEquationSystems();

        const std::filesystem::path trajectory_directory(
            input_db->getStringWithDefault("TRAJECTORY_DIRECTORY", "trajectory"));
        if (IBTK::IBTK_MPI::getRank() == 0) std::filesystem::create_directories(trajectory_directory);
        IBTK::IBTK_MPI::barrier();
        writeConnectivity(membrane_mesh, material_coordinates, trajectory_directory);

        int step = time_integrator->getIntegratorStep();
        double time = time_integrator->getIntegratorTime();
        double reference_area = 0.0;
        double reference_volume = 0.0;
        writeTrajectoryFrame(*equation_systems,
                             step,
                             time,
                             trajectory_directory,
                             reference_area,
                             reference_volume);

        if (dump_visualization)
        {
            time_integrator->setupPlotData();
            visit_writer->writePlotData(patch_hierarchy, step, time);
        }

        const int trajectory_interval = input_db->getIntegerWithDefault("TRAJECTORY_DUMP_INTERVAL", 5);
        while (!IBTK::rel_equal_eps(time, time_integrator->getEndTime()) && time_integrator->stepsRemaining())
        {
            const double dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            time += dt;
            ++step;
            if (step % trajectory_interval == 0 || !time_integrator->stepsRemaining())
                writeTrajectoryFrame(*equation_systems,
                                     step,
                                     time,
                                     trajectory_directory,
                                     reference_area,
                                     reference_volume);
            if (dump_visualization &&
                (step % visualization_interval == 0 || !time_integrator->stepsRemaining()))
            {
                time_integrator->setupPlotData();
                visit_writer->writePlotData(patch_hierarchy, step, time);
            }
        }

        for (auto* boundary_condition : velocity_boundary_conditions) delete boundary_condition;
    }
    return 0;
}
