/**
 * This file holds the implementation of the Sphere primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_SPHERE_H__
#define __DCA_SPHERE_H__

#include "Helpers.h"
#include "Primitives.h"

namespace DCA {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

double Sphere::compute_D(const primitive_t &other) const {
    auto d_sphere = [&](const Sphere &other) -> double {
        return (m_position - other.m_position).norm() -
               (m_radius + other.m_radius);
    };

    auto d_capsule = [&](const Capsule &other) -> double {
        VectorXd P(9);
        P << getPosition(), other.getStartPosition(), other.getEndPosition();
        return SphereCapsuleDistanceHelper::compute_D(P, getRadius(),
                                                      other.getRadius());
    };

    return std::visit(overloaded{d_sphere, d_capsule}, other);
}

void Sphere::compute_dDdP(VectorXd &dDdP, const primitive_t &other) const {
    auto dDdP_sphere = [&](const Sphere &other) -> void {
        Vector3d v = m_position - other.m_position;
        double v_norm = v.norm();

        if (v_norm < EPSILON) {
            v_norm = EPSILON;
        }

        dDdP = v / v_norm;
    };

    auto dDdP_capsule = [&](const Capsule &other) -> void {
        VectorXd P(9);
        P << getPosition(), other.getStartPosition(), other.getEndPosition();
        VectorXd dDdP_full;
        SphereCapsuleDistanceHelper::compute_dDdP(dDdP_full, P);
        dDdP = dDdP_full.head(3);
    };

    std::visit(overloaded{dDdP_sphere, dDdP_capsule}, other);
}

void Sphere::compute_d2DdP2(MatrixXd &d2DdP2, const primitive_t &other) const {
    auto d2DdP2_sphere = [&](const Sphere &other) -> void {
        Vector3d v = m_position - other.m_position;
        double v_norm = v.norm();

        if (v_norm < EPSILON) {
            v_norm = EPSILON;
        }
        d2DdP2.resize(3, 3);
        d2DdP2 =
            (Matrix3d::Identity() * v_norm - (v * v.transpose()) / v_norm) /
            (v_norm * v_norm);
    };

    auto d2DdP2_capsule = [&](const Capsule &other) -> void {
        VectorXd P(9);
        P << getPosition(), other.getStartPosition(), other.getEndPosition();
        MatrixXd d2DdP2_full;
        SphereCapsuleDistanceHelper::compute_d2DdP2(d2DdP2_full, P);
        d2DdP2 = d2DdP2_full.block(0, 0, 3, 3);
    };

    std::visit(overloaded{d2DdP2_sphere, d2DdP2_capsule}, other);
}
void Sphere::compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                                  const primitive_t &other) const {}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

}  // namespace DCA

#endif /* __DCA_SPHERE_H__ */