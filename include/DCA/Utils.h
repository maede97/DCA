/**
 * This file holds utils used all over this library.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_UTILS_H__
#define __DCA_UTILS_H__

#include <Eigen/Core>
#include <variant>
#include <vector>

namespace DCA {

// Forward definitions
class Sphere;
class Capsule;

// All possible primitves
using primitive_t = std::variant<Sphere, Capsule>;

// Easier access to a 3x3 matrix
using Eigen::Matrix3d;

// Easier access to a vector of length 3
using Eigen::Vector3d;

// Easier access to a vector of length 6
using Vector6d = Eigen::Matrix<double, 6, 1>;

// Easier access to any dynamic matrix
using Eigen::MatrixXd;

// Easier access to any dynamic vector
using Eigen::VectorXd;

// Easier access to a pair (corresponding of two indices)
using pair_t = std::pair<size_t, size_t>;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Helper for std::visit
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
#endif

/**
 * Simple helper to compute the sigmoid (with a scale and a shift)
 * @param x The value to take the sigmoid from.
 * @param scale The scale which is in the expontent.
 * @param shift The shift to apply to x before taking the sigmoid.
 */
inline double sigmoid(double x, double scale, double shift = 0.5) {
    return 1.0 / (1.0 + std::exp(-1.0 * scale * (x - shift)));
}

#define EPSILON 1e-8

/**
 * Helper for compute_d2DdP2_other
 */

/**
 * This class is the base class for all primitives.
 * with respect to ALL other primitives.
 */
class PrimitiveBase {
public:
    /**
     * Compute the distance to another object.
     * @param[in] other The other primitive to check against.
     * @return The distance to the other object.
     * @attention The primitive positions have to be in world coordinates!
     */
    virtual double compute_D(const primitive_t& other) const = 0;

    /**
     * Compute the first derivative (gradient)
     * of the distance function to the other object.
     * The gradient is taken with respect to this object,
     * with the other being non-movable.
     * @param[out] dDdP The gradient of the distance function.
     * @param[in] other The other primitive to check against.
     * @attention The primitive positions have to be in world coordinates!
     */
    virtual void compute_dDdP(VectorXd& dDdP, const primitive_t& other) const {
        /**
         * @todo
         */
    }

    /**
     * Compute the second derivative (hessian)
     * of the distance function to the other object.
     * The hessian is taken with respect to this object,
     * with the other being non-movable.
     * @param[out] d2DdP2 The hessian of the distance function.
     * @param[in] other The other primitive to check against.
     * @attention The primitive positions have to be in world coordinates!
     */
    virtual void compute_d2DdP2(MatrixXd& d2DdP2,
                                const primitive_t& other) const {
        /**
         * @todo
         */
    }

    /**
     * Compute the second derivative of the distance with respect to the others parameters.
     * @param[out] d2DdP2_other The hessian of the distance function (with respect to the other primitive)
     * @param[in] other The other primitive to check against.
     * @attention The primitive positions have to be in world coordinates!
     */
    virtual void compute_d2DdP2_other(MatrixXd& d2DdP2_other,
                                      const primitive_t& other) const {
        /**
         * @todo
         */
    }
};

/**
 * This class represents an objective used in Sensitivity Analysis.
 */
class SensitivityObjective {
    /**
     * Computes the distance for a given P and X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual double compute_D(const VectorXd& P, const VectorXd& X) const = 0;

    /**
     * Computes the first derivative of the distance with respect to X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dDdX(VectorXd& dDdX, const VectorXd& P,
                              const VectorXd& X) const = 0;

    /**
     * Computes the second derivative of the distance with respect to X.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdX2(MatrixXd& d2DdX2, const VectorXd& P,
                                const VectorXd& X) const = 0;

    /**
     * Computes the first derivative of the distance with respect to P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dDdP(VectorXd& dDdP, const VectorXd& P,
                              const VectorXd& X) const = 0;
    /**
     * Computes the second derivative of the distance with respect to P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdP2(MatrixXd& d2DdP2, const VectorXd& P,
                                const VectorXd& X) const = 0;
    /**
     * Computes the second derivative of the distance with respect to X and P.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2DdXdP(MatrixXd& d2DdXdP, const VectorXd& P,
                                 const VectorXd& X) const = 0;
};

}  // namespace DCA

#endif /* __DCA_UTILS_H__ */