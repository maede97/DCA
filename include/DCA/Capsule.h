/**
 * This file holds the implementation of the Capsule primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_CAPSULE_H__
#define __DCA_CAPSULE_H__

#include "Helpers.h"
#include "Primitives.h"

namespace DCA {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
double Capsule::compute_D(const primitive_t& other) const {
    auto d_capsule = [&](const Capsule& other) -> double {
        // Build P vector
        VectorXd P(12);
        P << getStartPosition(), getEndPosition(), other.getStartPosition(),
            other.getEndPosition();

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
        P << getStartPosition(), getEndPosition(), other.getStartPosition(),
            other.getEndPosition();

        CapsuleDistanceHelper::compute_dDdP(dDdP, P);
    };

    auto dDdP_sphere = [&](const Sphere& other) -> void {
        VectorXd P(9);
        P << other.getPosition(), getStartPosition(), getEndPosition();

        VectorXd dDdP_full;
        SphereCapsuleDistanceHelper::compute_dDdP(dDdP_full, P);
        dDdP = dDdP_full.tail(6);
    };

    std::visit(overloaded{dDdP_capsule, dDdP_sphere}, other);
}

void Capsule::compute_d2DdP2(MatrixXd& d2DdP2, const primitive_t& other) const {
    auto d2DdP2_capsule = [&](const Capsule& other) -> void {
        // Build the P vector
        VectorXd P(12);
        P << getStartPosition(), getEndPosition(), other.getStartPosition(),
            other.getEndPosition();

        CapsuleDistanceHelper::compute_d2DdP2(d2DdP2, P);
    };
    auto d2DdP2_sphere = [&](const Sphere& other) -> void {
        VectorXd P(9);
        P << other.getPosition(), getStartPosition(), getEndPosition();

        MatrixXd d2DdP2_full;
        SphereCapsuleDistanceHelper::compute_d2DdP2(d2DdP2_full, P);
        d2DdP2 = d2DdP2_full.block(3, 3, 6, 6);
    };

    std::visit(overloaded{d2DdP2_capsule, d2DdP2_sphere}, other);
}
void Capsule::compute_d2DdP2_other(MatrixXd& d2DdP2_other,
                                   const primitive_t& other) const {}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

}  // namespace DCA

#endif /* __DCA_CAPSULE_H__ */