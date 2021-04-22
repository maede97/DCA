#include <DCA/Newton.h>

#include <cmath>  // for std::isfinite

namespace DCA {

NewtonOptimizer::NewtonOptimizer(double solverResidual,
                                 unsigned int maxLineSearchIterations)
    : m_solverResidual(solverResidual),
      m_maxLineSearchIterations(maxLineSearchIterations) {}

bool NewtonOptimizer::optimize(NewtonObjective& objective, const VectorXd& P,
                               VectorXd& x, unsigned int maxIterations) {
    {
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
}

void NewtonOptimizer::computeSearchDirection(const VectorXd& P,
                                             NewtonObjective& objective) {
    objective.compute_dOdX(m_gradient, P, m_x_tmp);
    objective.compute_d2OdX2(m_hessian, P, m_x_tmp);

    solveLinearSystem(m_searchDir, m_hessian, m_gradient);

    if (m_useDynamicRegularization)
        applyDynamicRegularization(m_searchDir, m_hessian, m_gradient);
}

bool NewtonOptimizer::doLineSearch(const VectorXd& P,
                                   NewtonObjective& objective) {
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

void NewtonOptimizer::solveLinearSystem(VectorXd& y, const MatrixXd& A,
                                        const VectorXd& x) {
    y = A.colPivHouseholderQr().solve(x);
}

void NewtonOptimizer::applyDynamicRegularization(VectorXd& y, MatrixXd& A,
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

}  // namespace DCA