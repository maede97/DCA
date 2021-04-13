#ifndef __DCA_NEWTON_H__
#define __DCA_NEWTON_H__

#include <cmath>  // for std::isfinite
#include <Eigen/Dense>
#include <iostream>

#include "utils.h"

namespace DCA {

/**
 * @todo: 
 * don't store class members, but
 * rather give by reference/value
 */
class NewtonOptimizer {
public:
    NewtonOptimizer(double solverResidual = 1e-5,
                    unsigned int maxLineSearchIterations = 15)
        : m_solverResidual(solverResidual),
          m_maxLineSearchIterations(maxLineSearchIterations) {}

    ~NewtonOptimizer() {}

    // also add P vector here
    bool optimize(Objective& objective, const VectorXd& P, VectorXd& x,
                  unsigned int maxIterations = 100) {
        m_x_tmp = x;
        bool betterSolutionFound = false;
        bool converged = false;

        for (int i = 0; i < maxIterations; i++) {
            objective.preOptimizationStep(P, m_x_tmp);

            computeSearchDirection(P, objective);
            betterSolutionFound = doLineSearch(P, objective);

            objective.postOptimizationStep(P, m_x_tmp);

            if (m_gradient.norm() < m_solverResidual) {
                converged = true;
                break;
            }
        }

        if (betterSolutionFound) x = m_x_tmp;
        return converged;
    }

private:
    void computeSearchDirection(const VectorXd& P, Objective& objective) {
        objective.compute_dOdX(m_gradient, P, m_x_tmp);
        objective.compute_d2OdX2(m_hessian, P, m_x_tmp);

        solveLinearSystem(m_searchDir, m_hessian, m_gradient);

        if (m_useDynamicRegularization)
            applyDynamicRegularization(m_searchDir, m_hessian, m_gradient);
    }
    bool doLineSearch(const VectorXd& P, Objective& objective) {
        if (m_maxLineSearchIterations < 1) {
            m_x_tmp = m_x_tmp - m_searchDir * m_lineSearchStartValue;
            return true;
        }

        double alpha = m_lineSearchStartValue;
        VectorXd xc(m_x_tmp);
        double initialObjectiveValue = objective.compute_O(P, xc);

        for (int j = 0; j < m_maxLineSearchIterations; j++) {
            m_x_tmp = xc - m_searchDir * alpha;
            objective.postLineSearchStep(P, m_x_tmp);
            m_objectiveValue = objective.compute_O(P, m_x_tmp);

            if (!std::isfinite(m_objectiveValue))
                m_objectiveValue = initialObjectiveValue + 1.;

            if (m_objectiveValue > initialObjectiveValue)
                alpha *= 0.5;
            else
                return true;
        }
        return false;
    }

    static void solveLinearSystem(VectorXd& y, const MatrixXd& A,
                                  const VectorXd& b) {
        y = A.colPivHouseholderQr().solve(b);
    }

    static void applyDynamicRegularization(VectorXd& y, MatrixXd& A,
                                           const VectorXd& x) {
        double dotProduct = y.dot(x);
        if (dotProduct <= 0. && x.squaredNorm() > 0) {
            VectorXd stabRegularizer(A.rows());
            stabRegularizer.setZero();

            double currStabValue = 1e-4;
            for (int i = 0; i < 10; i++) {
                stabRegularizer.setConstant(currStabValue);
                A += stabRegularizer.asDiagonal();
                currStabValue *= 10.;

                solveLinearSystem(y, A, x);

                dotProduct = y.dot(x);
                if (dotProduct > 0.) break;
            }
        }
    }

private:
    double m_solverResidual;
    unsigned int m_maxLineSearchIterations;
    double m_lineSearchStartValue = 1.0;
    bool m_useDynamicRegularization = true;

    double m_objectiveValue;
    VectorXd m_x_tmp, m_searchDir, m_gradient;
    MatrixXd m_hessian;
};

}  // namespace DCA

#endif /* __DCA_NEWTON_H__ */