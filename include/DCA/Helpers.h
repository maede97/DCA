/**
 * This file holds helpers to compute primitive distances and derivatives.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_HELPERS_H__
#define __DCA_HELPERS_H__

#include "Newton.h"
#include "SoftUpperLimitConstraint.h"
#include "Utils.h"

namespace DCA {

/**
 * This objective is used to solve the capsule-distance problem:
 * A capsule is represented as a line between two points.
 * We search now the shortest distance between two capsules, i.e. two lines.
 * We solve it using Newont's Method (hence the inheritance of Objective).
 * Furthermore, we provide all derivatives for Sensitivity Analysis.
 */
class CapsuleDistanceObjective : public NewtonObjective,
                                 public SensitivityObjective {
public:
    /**
     * @copydoc NewtonObjective::compute_O
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    double compute_O(const VectorXd& P, const VectorXd& X) const override;

    /**
     * @copydoc NewtonObjective::compute_dOdX
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_dOdX(VectorXd& dOdX, const VectorXd& P,
                      const VectorXd& X) const override;

    /**
     * @copydoc NewtonObjective::compute_d2OdX2
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2OdX2(MatrixXd& d2OdX2, const VectorXd& P,
                        const VectorXd& X) const override;

    /**
     * @copydoc SensitivityObjective::compute_dDdX
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_dDdX(VectorXd& dOdX, const VectorXd& P,
                      const VectorXd& X) const override;
    /**
     * @copydoc SensitivityObjective::compute_d2DdX2
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2DdX2(MatrixXd& d2OdX2, const VectorXd& P,
                        const VectorXd& X) const override;

    /**
     * @copydoc SensitivityObjective::compute_D
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    double compute_D(const VectorXd& P, const VectorXd& X) const override;

    /**
     * @copydoc SensitivityObjective::compute_dDdP
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_dDdP(VectorXd& dDdP, const VectorXd& P,
                      const VectorXd& X) const override;

    /**
     * @copydoc SensitivityObjective::compute_d2DdP2
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2DdP2(MatrixXd& d2DdP2, const VectorXd& P,
                        const VectorXd& X) const override;

    /**
     * @copydoc SensitivityObjective::compute_d2DdXdP
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2DdXdP(MatrixXd& d2DdXdP, const VectorXd& P,
                         const VectorXd& X) const override;

private:
    double regularizerWeight = 0.1;  ///< weight for the regularizer
    double constraintWeight = 10.0;  ///< weight for the constraint

    /// 2: Both are regularized / 1: only one is regularized
    int maxRegularizerIndex = 2;

    /// The constraint with constrains the value
    SoftUpperLimitConstraint sulc = SoftUpperLimitConstraint(0.0, 1.0, 0.001);
};

/**
 * This class helps you compute collision avoidance for capsule vs. capsule.
 */
class CapsuleDistanceHelper {
public:
    /**
     * Computes the first derivative of X with respect to P.
     */
    static void compute_dXdP(MatrixXd& dXdP, const VectorXd& P,
                             const VectorXd& X);

    /**
     * Computes the distance between both capsules.
     * @param P Stacked end points, where [0 -> 2] and [3 -> 5] are from the first capsule
     * and [6 -> 8] and [9 -> 11] are from the second capsule.
     * @param r12 The radius of the first capsule
     * @param r34 The radius of the second capsule
     * @return The shortest distance between both capsules.
     */
    static double compute_D(const VectorXd& P, const double& r12,
                            const double& r34);

    /**
     * Computes the second derivative of the distance with respect to
     * the points of the first capsule.
     */
    static void compute_dDdP(VectorXd& dDdP, const VectorXd& P);

    /**
     * Computes the first derivative of the distance with respect to
     * the points of the first capsule.
     */
    static void compute_d2DdP2(MatrixXd& d2DdP2, const VectorXd& P);

private:
    /**
     * Compute X based on the P variable
     */
    static void solveForX(const VectorXd& P, VectorXd& X);

    /**
     * Compute the closest points on two lines based on P.
     */
    static void computeClosestPointOnLines(Vector3d& P12, Vector3d& P34,
                                           const VectorXd& P);

    /**
     * The objective to use
     */
    static CapsuleDistanceObjective objective;
};

/**
 * This helper is used to compute sphere-capsule distances.
 */
class SphereCapsuleDistanceHelper {
public:
    /**
     * Compute the distance between a sphere and a capsule.
     * @param P Is stacked with one Vector3d (Sphere Position) and two Vector3ds (Capsule end points)
     * @param rSphere The radius of the sphere.
     * @param rCapsule The radius of the capsule.
     */
    static double compute_D(const VectorXd& P, const double& rSphere,
                            const double& rCapsule);

    /**
     * Compute the first derivative distance with respect to the sphere AND capsule.
     * @param P Is stacked with one Vector3d (Sphere Position) and two Vector3ds (Capsule end points)
     * @attention Computes the derivative with respect to the sphere AND capsule.
     */
    static void compute_dDdP(VectorXd& dDdP, const VectorXd& P);

    /**
     * Compute the second derivative of the distance with respect to the sphere AND capsule
     * @param P Is stacked with one Vector3d (Sphere Position) and two Vector3ds (Capsule end points)
     * @attention Computes the derivative with respect to the sphere AND capsule.
     */
    static void compute_d2DdP2(MatrixXd& d2DdP2, const VectorXd& P);

private:
    static Vector3d computeClosestPointOnLine(const VectorXd& P);

    static double compute_t(const VectorXd& P);

    static double sigScale();
};

}  // namespace DCA

#endif /* __DCA_HELPERS_H__ */