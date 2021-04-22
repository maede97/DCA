/**
 * This file holds a simple newton optimizer.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_NEWTON_H__
#define __DCA_NEWTON_H__

#include <Eigen/Dense>

#include "Utils.h"

namespace DCA {

/**
 * This class represents an objective, which can be used together with the NewtonMinimizer.
 */
class NewtonObjective {
public:
    /**
     * Computes the Objective value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual double compute_O(const VectorXd& P, const VectorXd& X) const = 0;
    /**
     * Computes the derivative of the objective value with respect to X value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_dOdX(VectorXd& dOdX, const VectorXd& P,
                              const VectorXd& X) const = 0;
    /**
     * Computes the second derivative of the objective value with respect to X value of this.
     * @param P Some parameters for the objective.
     * @param X The solving variable.
     */
    virtual void compute_d2OdX2(MatrixXd& d2OdX2, const VectorXd& P,
                                const VectorXd& X) const = 0;

    virtual void preOptimizationStep(const VectorXd& P, const VectorXd& X) {}
    virtual void postOptimizationStep(const VectorXd& P, const VectorXd& X) {}
    virtual void postLineSearchStep(const VectorXd& P, const VectorXd& X) {}

public:
    double weight;
};

/**
 * @todo: 
 * don't store class members, but
 * rather give by reference/value
 */
class NewtonOptimizer {
public:
    NewtonOptimizer(double solverResidual = 1e-5,
                    unsigned int maxLineSearchIterations = 15);

    ~NewtonOptimizer() {}

    /**
     * Optimize the objective.
     * @return true if the solver converged, false otherwise.
     */
    bool optimize(NewtonObjective& objective, const VectorXd& P, VectorXd& x,
                  unsigned int maxIterations = 100);

private:
    /**
     * Compute the search direction by searching for the gradient direction
     */
    void computeSearchDirection(const VectorXd& P, NewtonObjective& objective);

    /**
     * Perform line search on the objective
     */
    bool doLineSearch(const VectorXd& P, NewtonObjective& objective);

    /**
     * Solve the system A * y = x for y.
     */
    static void solveLinearSystem(VectorXd& y, const MatrixXd& A,
                                  const VectorXd& x);

    /**
     * Apply dynamic regularization on the linear system A * y = x
     */
    static void applyDynamicRegularization(VectorXd& y, MatrixXd& A,
                                           const VectorXd& x);

private:
    ///< residual of the solver
    double m_solverResidual;
    ///< how many line search steps should be done
    unsigned int m_maxLineSearchIterations;
    ///< the starting value of the line search
    double m_lineSearchStartValue = 1.0;
    ///< whether to use dynamic regularization
    bool m_useDynamicRegularization = true;

    ///< the current objective value
    double m_objectiveValue;
    ///< some private members, could also be given via parameters
    VectorXd m_x_tmp, m_searchDir, m_gradient;
    ///< and the hessian
    MatrixXd m_hessian;
};

}  // namespace DCA

#endif /* __DCA_NEWTON_H__ */