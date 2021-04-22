/**
 * This file holds the implementation for each shape primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_PRIMITIVES_H__
#define __DCA_PRIMITIVES_H__

#include "Utils.h"

namespace DCA {

/**
 * Shape Primitive: Sphere
 */
class Sphere : public PrimitiveBase {
public:
    Sphere(const Vector3d &position, const double &radius) : m_radius(radius) {
        m_parameters.resize(3);
        m_parameters << position;
    }

    Sphere(const VectorXd params) : PrimitiveBase(params) {}

    double compute_D(const primitive_t &other) const override;
    void compute_dDdP(VectorXd &dDdP, const primitive_t &other) const override;
    void compute_d2DdP2(MatrixXd &d2DdP2,
                        const primitive_t &other) const override;
    void compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                              const primitive_t &other) const override;

    // Helpers
    Vector3d getPosition() const { return m_parameters.head(3); }
    double getRadius() const { return m_radius; }

private:
    double m_radius;
};

class Capsule : public PrimitiveBase {
public:
    Capsule(const Vector3d &startPosition, const Vector3d &endPosition,
            const double &radius) : m_radius(radius) {
        m_parameters.resize(6);
        m_parameters << startPosition, endPosition;
    }

    Capsule(const VectorXd params) : PrimitiveBase(params) {}

    double compute_D(const primitive_t &other) const override;
    void compute_dDdP(VectorXd &dDdP, const primitive_t &other) const;
    void compute_d2DdP2(MatrixXd &d2DdP2,
                        const primitive_t &other) const override;
    void compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                              const primitive_t &other) const override;

    // Helpers
    Vector3d getStartPosition() const { return m_parameters.head(3); }
    Vector3d getEndPosition() const { return m_parameters.tail(3); }
    double getRadius() const { return m_radius; }

private:
    double m_radius;
};

// ========== CUSTOM FUNCTIONS ==========
// The following three functions are for easier access from outside
// such that the user does not need to write custom std::visit calls.

/**
 * Compute the distance of a given pair (indices) from the primitives
 */
double compute_D(const pair_t &pair,
                 const std::vector<primitive_t> &primitives) {
    return std::visit([](const auto &cp1,
                         const primitive_t &cp2) { return cp1.compute_D(cp2); },
                      primitives.at(pair.first), primitives.at(pair.second));
}

/**
 * Computes the gradient of a given pair and primitives
 */
void compute_dDdP(Eigen::VectorXd &dDdP, const pair_t &pair,
                  const std::vector<primitive_t> &primitives) {
    std::visit([&dDdP](const auto &cp1,
                       const primitive_t &cp2) { cp1.compute_dDdP(dDdP, cp2); },
               primitives.at(pair.first), primitives.at(pair.second));
}

/**
 * Computes the hessian of a given pair and primitives
 */
void compute_d2DdP2(Eigen::MatrixXd &d2DdP2, const pair_t &pair,
                    const std::vector<primitive_t> &primitives) {
    std::visit(
        [&d2DdP2](const auto &cp1, const primitive_t &cp2) {
            cp1.compute_d2DdP2(d2DdP2, cp2);
        },
        primitives.at(pair.first), primitives.at(pair.second));
}

}  // namespace DCA
#endif /* __DCA_PRIMITIVES_H__ */