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

/**
 * Shape Primitive: Sphere
 */
class Sphere : public PrimitiveBase {
public:
    Sphere(const Vector3d &position, const double &radius);

    Sphere(const VectorXd params) : PrimitiveBase(params) {}

    double compute_D(const primitive_t &other) const override;
    void compute_dDdP(VectorXd &dDdP, const primitive_t &other) const override;
    void compute_d2DdP2(MatrixXd &d2DdP2,
                        const primitive_t &other) const override;
    void compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                              const primitive_t &other) const override;

    // Helpers
    Vector3d getPosition() const { return m_parameter.head(3); }
    double getRadius() const { return m_radius; }

private:
    double m_radius;
};

}  // namespace DCA

#endif /* __DCA_SPHERE_H__ */
