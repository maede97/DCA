#ifndef __DCA_NEWTON_H__
#define __DCA_NEWTON_H__

#include <cmath>  // for std::isfinite
#include <eigen3/Eigen/SparseLU>
#include <iostream>
#include "utils.h"

namespace DCA {

class Objective {
public:
    virtual double evaluate(const VectorXd& x) const = 0;
    virtual VectorXd gradient(const VectorXd& x) const = 0;
    virtual SparseMatrixd hessian(const VectorXd& x) const = 0;

    virtual void preOptimizationStep(const VectorXd& x) {}
    virtual void postOptimizationStep(const VectorXd& x) {}
    virtual void postLineSearchStep(const VectorXd& x) {}

public:
    double weight;
};

class NewtonOptimizer {
public:
    NewtonOptimizer(double solverResidual = 1e-5,
                    unsigned int maxLineSearchIterations = 15)
        : m_solverResidual(solverResidual),
          m_maxLineSearchIterations(maxLineSearchIterations) {}

    ~NewtonOptimizer() {}

    bool optimize(Objective* objective, VectorXd& x,
                  unsigned int maxIterations = 100) {
        m_x_tmp = x;
        bool betterSolutionFound = false;
        bool converged = false;

        for (int i = 0; i < maxIterations; i++) {
            objective->preOptimizationStep(m_x_tmp);

            computeSearchDirection(objective);
            betterSolutionFound = doLineSearch(objective);

            objective->postOptimizationStep(m_x_tmp);

            if (m_gradient.norm() < m_solverResidual) {
                converged = true;
                break;
            }
        }

        if (betterSolutionFound) x = m_x_tmp;
        return converged;
    }

private:
    void computeSearchDirection(Objective* objective) {
        m_gradient = objective->gradient(m_x_tmp);
        m_hessian = objective->hessian(m_x_tmp);

        solveLinearSystem(m_searchDir, m_hessian, m_gradient);

        if (m_useDynamicRegularization)
            applyDynamicRegularization(m_searchDir, m_hessian, m_gradient);
    }
    bool doLineSearch(Objective* objective) {
        if (m_maxLineSearchIterations < 1) {
            m_x_tmp = m_x_tmp - m_searchDir * m_lineSearchStartValue;
            return true;
        }

        double alpha = m_lineSearchStartValue;
        VectorXd xc(m_x_tmp);
        double initialObjectiveValue = objective->evaluate(xc);

        for (int j = 0; j < m_maxLineSearchIterations; j++) {
            std::cout << m_x_tmp.size() << std::endl;
            std::cout << m_searchDir.size() << std::endl;
            m_x_tmp = xc - m_searchDir * alpha;
            objective->postLineSearchStep(m_x_tmp);
            m_objectiveValue = objective->evaluate(m_x_tmp);

            if (!std::isfinite(m_objectiveValue))
                m_objectiveValue = initialObjectiveValue + 1.;

            if (m_objectiveValue > initialObjectiveValue)
                alpha *= 0.5;
            else
                return true;
        }
        return false;
    }

    static void solveLinearSystem(VectorXd& y, const SparseMatrixd& A,
                                  const VectorXd& b) {
        Eigen::SimplicialLDLT<SparseMatrixd, Eigen::Lower> solver;
        solver.compute(A);
        y = solver.solve(b);
    }
    static void applyDynamicRegularization(VectorXd& y, SparseMatrixd& A,
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
    SparseMatrixd m_hessian;
};

}  // namespace DCA

#endif /* __DCA_NEWTON_H__ */