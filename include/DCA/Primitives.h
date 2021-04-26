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

// ========== CUSTOM FUNCTIONS ==========
// The following functions are for easier access from outside
// such that the user does not need to write custom std::visit calls.

/**
 * Compute the distance of a given pair (indices) from the primitives
 */
double compute_D(const pair_t &pair,
                 const std::vector<primitive_t> &primitives);

/**
 * Computes the gradient of a given pair and primitives
 */
void compute_dDdP(VectorXd &dDdP, const pair_t &pair,
                  const std::vector<primitive_t> &primitives);

/**
 * Computes the hessian of a given pair and primitives.
 */
void compute_d2DdP2(MatrixXd &d2DdP2, const pair_t &pair,
                    const std::vector<primitive_t> &primitives);

/**
 * Computes the other hessian of a given pair and primitives.
 */
void compute_d2DdP2_other(MatrixXd& d2DdP2_other, const pair_t& pair, const std::vector<primitive_t>& primitives);

/**
 * Get's a vector with all parameters that define a primitive (and that change the derivatives).
 * E.g. for a Sphere, this is only the position. For a Capsule, it's both the start and end points, stacked.
 */
void getParameter(VectorXd &parameter, const primitive_t primitive);

}  // namespace DCA
#endif /* __DCA_PRIMITIVES_H__ */