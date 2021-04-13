/**
 * This file holds the implementation for each shape primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_PRIMITIVES_H__
#define __DCA_PRIMITIVES_H__

#include "utils.h"

namespace DCA {

/**
 * Shape Primitive: Sphere
 */
class Sphere : public PrimitiveBase {
public:
    Sphere(const Vector3d &position, const double &radius)
        : m_position(position), m_radius(radius) {}

    double compute_D(const primitive_t &other) const override;
    void compute_dDdP(VectorXd &dDdP, const primitive_t &other) const override;
    void compute_d2DdP2(MatrixXd &d2DdP2,
                        const primitive_t &other) const override;
    void compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                              const primitive_t &other) const override;

    // Helpers
    Vector3d getPosition() const { return m_position; }
    double getRadius() const { return m_radius; }

private:
    Vector3d m_position;
    double m_radius;
};

class Capsule : public PrimitiveBase {
public:
    Capsule(const Vector3d &startPosition, const Vector3d &endPosition,
            const double &radius)
        : m_startPosition(startPosition),
          m_endPosition(endPosition),
          m_radius(radius) {}

    double compute_D(const primitive_t &other) const override;
    void compute_dDdP(VectorXd &dDdP, const primitive_t &other) const;
    void compute_d2DdP2(MatrixXd &d2DdP2,
                        const primitive_t &other) const override;
    void compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                              const primitive_t &other) const override;

    // Helpers
    Vector3d getStartPosition() const { return m_startPosition; }
    Vector3d getEndPosition() const { return m_endPosition; }
    double getRadius() const { return m_radius; }

private:
    Vector3d m_startPosition;
    Vector3d m_endPosition;
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