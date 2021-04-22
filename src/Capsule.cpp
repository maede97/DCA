#include <DCA/Capsule.h>
#include <DCA/Sphere.h>

namespace DCA {

Capsule::Capsule(const Vector3d& startPosition, const Vector3d& endPosition,
                 const double& radius)
    : m_radius(radius) {
    m_parameter.resize(6);
    m_parameter << startPosition, endPosition;
}

double Capsule::compute_D(const primitive_t& other) const {
    auto d_capsule = [&](const Capsule& other) -> double {
        // Build P vector
        VectorXd P(12);
        P << other.getStartPosition(), other.getEndPosition(),
            getStartPosition(), getEndPosition();

        // Use the CapsuleDistanceHelper to compute the distance
        return CapsuleDistanceHelper::compute_D(P, getRadius(),
                                                other.getRadius());
    };

    auto d_sphere = [&](const Sphere& other) -> double {
        VectorXd P(9);
        P << other.getPosition(), getStartPosition(), getEndPosition();
        return SphereCapsuleDistanceHelper::compute_D(P, getRadius(),
                                                      other.getRadius());
    };

    return std::visit(overloaded{d_capsule, d_sphere}, other);
}

void Capsule::compute_dDdP(VectorXd& dDdP, const primitive_t& other) const {
    auto dDdP_capsule = [&](const Capsule& other) -> void {
        // Build the P vector
        VectorXd P(12);
        P << other.getStartPosition(), other.getEndPosition(),
            getStartPosition(), getEndPosition();

        CapsuleDistanceHelper::compute_dDdP(dDdP, P);
    };

    auto dDdP_sphere = [&](const Sphere& other) -> void {
        VectorXd P(9);
        P << other.getPosition(), getStartPosition(), getEndPosition();

        VectorXd dDdP_full;
        SphereCapsuleDistanceHelper::compute_dDdP(dDdP_full, P);
        dDdP = dDdP_full.head(3);
    };

    std::visit(overloaded{dDdP_capsule, dDdP_sphere}, other);
}

void Capsule::compute_d2DdP2(MatrixXd& d2DdP2, const primitive_t& other) const {
    auto d2DdP2_capsule = [&](const Capsule& other) -> void {
        // Build the P vector
        VectorXd P(12);
        P << other.getStartPosition(), other.getEndPosition(),
            getStartPosition(), getEndPosition();

        CapsuleDistanceHelper::compute_d2DdP2(d2DdP2, P);
    };
    auto d2DdP2_sphere = [&](const Sphere& other) -> void {
        VectorXd P(9);
        P << other.getPosition(), getStartPosition(), getEndPosition();

        MatrixXd d2DdP2_full;
        SphereCapsuleDistanceHelper::compute_d2DdP2(d2DdP2_full, P);
        d2DdP2 = d2DdP2_full.block(0, 0, 3, 3);
    };

    std::visit(overloaded{d2DdP2_capsule, d2DdP2_sphere}, other);
}
void Capsule::compute_d2DdP2_other(MatrixXd& d2DdP2_other,
                                   const primitive_t& other) const {}

}  // namespace DCA