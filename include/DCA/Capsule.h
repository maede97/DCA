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

class Capsule : public PrimitiveBase {
public:
    Capsule(const Vector3d &startPosition, const Vector3d &endPosition,
            const double &radius);

    Capsule(const VectorXd params) : PrimitiveBase(params) {}

    double compute_D(const primitive_t &other) const override;
    void compute_dDdP(VectorXd &dDdP, const primitive_t &other) const;
    void compute_d2DdP2(MatrixXd &d2DdP2,
                        const primitive_t &other) const override;
    void compute_d2DdP2_other(MatrixXd &d2DdP2_other,
                              const primitive_t &other) const override;

    // Helpers
    Vector3d getStartPosition() const { return m_parameter.head(3); }
    Vector3d getEndPosition() const { return m_parameter.tail(3); }
    double getRadius() const { return m_radius; }

private:
    double m_radius;
};

}  // namespace DCA

#endif /* __DCA_CAPSULE_H__ */