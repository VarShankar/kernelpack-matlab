#ifndef KERNELPACK_RBC_SURFACE_METHOD_H
#define KERNELPACK_RBC_SURFACE_METHOD_H

#include <ibamr/IBFESurfaceMethod.h>

#include <array>
#include <map>
#include <string>
#include <vector>

namespace KernelPackRBC
{
struct MembraneParameters
{
    double shear_modulus = 0.5;
    double bending_modulus = 2.5e-4;
    double area_modulus = 5.0;
    double volume_modulus = 20.0;
};

class RBCSurfaceMethod : public IBAMR::IBFESurfaceMethod
{
public:
    using IBAMR::IBFESurfaceMethod::IBFESurfaceMethod;

    static const std::string MEMBRANE_FORCE_SYSTEM_NAME;

    void initializeMembraneModel(const MembraneParameters& parameters, unsigned int part = 0);

    void computeLagrangianForce(double data_time) override;

private:
    struct Face
    {
        std::array<libMesh::dof_id_type, 3> node_ids;
    };

    struct Edge
    {
        libMesh::dof_id_type i;
        libMesh::dof_id_type j;
        libMesh::dof_id_type k;
        libMesh::dof_id_type l;
        double reference_length;
        double reference_dihedral;
    };

    void buildReferenceGeometry(unsigned int part);
    void updateMembraneForce(double data_time, unsigned int part);

    static double dihedralAngle(const std::array<libMesh::Point, 4>& x);

    MembraneParameters d_parameters;
    std::vector<Face> d_faces;
    std::vector<Edge> d_edges;
    std::vector<libMesh::Point> d_reference_positions;
    std::vector<double> d_reference_lumped_areas;
    double d_reference_area = 0.0;
    double d_reference_volume = 0.0;
    unsigned int d_part = 0;
    bool d_membrane_initialized = false;
};

void interpolateMembraneForce(
    libMesh::VectorValue<double>& force,
    const libMesh::VectorValue<double>& current_normal,
    const libMesh::VectorValue<double>& reference_normal,
    const libMesh::TensorValue<double>& deformation_gradient,
    const libMesh::Point& current_position,
    const libMesh::Point& reference_position,
    libMesh::Elem* element,
    unsigned short side,
    const std::vector<const std::vector<double>*>& variable_data,
    const std::vector<const std::vector<libMesh::VectorValue<double>>*>& variable_gradient_data,
    double data_time,
    void* context);
} // namespace KernelPackRBC

#endif
