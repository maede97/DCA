/**
 * This file holds a soft upper limit constraint.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_SOFTUPPERLIMITCONSTRAINT_H__
#define __DCA_SOFTUPPERLIMITCONSTRAINT_H__

namespace DCA {

class SoftUpperLimitConstraint {
private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon;
    double limit = 0;

public:
    SoftUpperLimitConstraint(double l, double stiffness, double epsilon);

    virtual ~SoftUpperLimitConstraint() {}

    void setLimit(double l) { limit = l; }
    void setEpsilon(double eps) { epsilon = eps; }

    double compute_F(double x) const;

    double compute_dFdX(double x) const;

    double compute_d2FdX2(double x) const;
};

}  // namespace DCA

#endif /* __DCA_SOFTUPPERLIMITCONSTRAINT_H__ */