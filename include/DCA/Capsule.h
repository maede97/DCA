/**
 * This file holds the implementation of the Capsule primitive.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_CAPSULE_H__
#define __DCA_CAPSULE_H__

#include "newton.h"
#include "primitives.h"

namespace DCA {

/**
 * This objective is used to solve the capsule-distance problem:
 * A capsule is represented as a line between two points.
 * We search now the shortest distance between two capsules, i.e. two lines.
 * We solve it using Newont's Method (hence the inheritance of Objective).
 * Furthermore, we provide all derivatives for Sensitivity Analysis.
 */
class CapsuleDistanceObjective : public Objective, public SensitivityObjective {
public:
    /**
     * @copydoc Objective::compute_O
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    double compute_O(const VectorXd& P, const VectorXd& X) const override {
        Vector3d P12 = P.head(3) + X[0] * (P.segment(3, 3) - P.head(3));
        Vector3d P34 = P.segment(6, 3) + X[1] * (P.tail(3) - P.segment(6, 3));
        return (P12 - P34).norm();
    }
    /**
     * @copydoc Objective::compute_dOdX
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_dOdX(VectorXd& dOdX, const VectorXd& P,
                      const VectorXd& X) const override {
        Vector3d P1 = P.segment(0, 3);
        Vector3d P2 = P.segment(3, 3);
        Vector3d P3 = P.segment(6, 3);
        Vector3d P4 = P.segment(9, 3);

        Vector3d P12 = P1 + X[0] * (P2 - P1);
        Vector3d P34 = P3 + X[1] * (P4 - P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;

        dOdX.resize(2);
        dOdX[0] = (P2 - P1).transpose() * (v / v_norm);
        dOdX[1] = -(P4 - P3).transpose() * (v / v_norm);
    }
    /**
     * @copydoc Objective::compute_d2OdX2
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2OdX2(MatrixXd& d2OdX2, const VectorXd& P,
                        const VectorXd& X) const override {
        Vector3d P1 = P.segment(0, 3);
        Vector3d P2 = P.segment(3, 3);
        Vector3d P3 = P.segment(6, 3);
        Vector3d P4 = P.segment(9, 3);

        Vector3d P12 = P1 + X[0] * (P2 - P1);
        Vector3d P34 = P3 + X[1] * (P4 - P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;
        d2OdX2.resize(2, 2);
        d2OdX2(0, 0) = (double)((P2 - P1).transpose() *
                                ((P2 - P1) * v_norm -
                                 (v * (P2 - P1).transpose() * v / v_norm))) /
                       (v_norm * v_norm);
        d2OdX2(1, 0) = (double)((P4 - P3).transpose() *
                                ((P2 - P1) * v_norm -
                                 (v * (P2 - P1).transpose() * v / v_norm))) /
                       (v_norm * v_norm) * -1.0;
        d2OdX2(0, 1) = (double)((P2 - P1).transpose() *
                                ((P4 - P3) * v_norm -
                                 (v * (P4 - P3).transpose() * v / v_norm))) /
                       (v_norm * v_norm) * -1.0;
        d2OdX2(1, 1) = (double)((P4 - P3).transpose() *
                                ((P4 - P3) * v_norm -
                                 (v * (P4 - P3).transpose() * v / v_norm))) /
                       (v_norm * v_norm);
    }

    /**
     * @copydoc SensitivityObjective::compute_D
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    double compute_D(const VectorXd& P, const VectorXd& X) const override {
        /**
         * @todo: Is this the one with regularizers? And then, what about getGradient
         */
        return compute_O(P, X);
    }

    /**
     * @copydoc SensitivityObjective::compute_dDdP
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_dDdP(VectorXd& dDdP, const VectorXd& P,
                      const VectorXd& X) const override {
        Vector3d P1 = P.segment(0, 3);
        Vector3d P2 = P.segment(3, 3);
        Vector3d P3 = P.segment(6, 3);
        Vector3d P4 = P.segment(9, 3);

        Vector3d P12 = P1 + X[0] * (P2 - P1);
        Vector3d P34 = P3 + X[1] * (P4 - P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;

        dDdP.resize(6);
        dDdP.segment(0, 3) = (1.0 - X[0]) * (v / v_norm);
        dDdP.segment(3, 3) = X[0] * (v / v_norm);
    }

    /**
     * @copydoc SensitivityObjective::compute_d2DdP2
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2DdP2(MatrixXd& d2DdP2, const VectorXd& P,
                        const VectorXd& X) const override {
        const double t12 = X[0];
        const double t34 = X[1];

        Vector3d P1 = P.segment(0, 3);
        Vector3d P2 = P.segment(3, 3);
        Vector3d P3 = P.segment(6, 3);
        Vector3d P4 = P.segment(9, 3);

        Vector3d P12 = P1 + t12 * (P2 - P1);
        Vector3d P34 = P3 + t34 * (P4 - P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;

        const Matrix3d I_vn = v_norm * Matrix3d::Identity();
        const Matrix3d v_T_v = v * (v / v_norm).transpose();
        const double v_norm_sq = v_norm * v_norm;

        double du_p, dv_p;
        Matrix3d ppD;
        double du = 1.0 - t12;

        d2DdP2.resize(6, 6);
        ////// pDpP1 = (1.0 - t12) * (v / v_norm)
        {  //--- dP1
            du_p = (1.0 - t12) * (1.0 - t12);
            dv_p = (1.0 - t12);
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            d2DdP2.block(0, 0, 3, 3) = ppD;
        }
        {  //--- dP2
            du_p = t12 * (1.0 - t12);
            dv_p = t12;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            d2DdP2.block(3, 0, 3, 3) = ppD;
            d2DdP2.block(0, 3, 3, 3) = ppD;
        }

        ////// pDpP2 = t12 * (v / v_norm)
        du = t12;
        {  //--- dP2
            du_p = t12 * t12;
            dv_p = t12;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            d2DdP2.block(3, 3, 3, 3) = ppD;
        }
    }

    /**
     * @copydoc SensitivityObjective::compute_d2DdXdP
     * @param P The four points, stacked.
     * @param X The current point on the two lines.
     */
    void compute_d2DdXdP(MatrixXd& d2DdXdP, const VectorXd& P,
                         const VectorXd& X) const override {
        Vector3d P1 = P.segment(0, 3);
        Vector3d P2 = P.segment(3, 3);
        Vector3d P3 = P.segment(6, 3);
        Vector3d P4 = P.segment(9, 3);

        Vector3d P12 = P1 + X[0] * (P2 - P1);
        Vector3d P34 = P3 + X[1] * (P4 - P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;
        Matrix3d I = Matrix3d::Identity();

        d2DdXdP.resize(2, 6);

        auto setEntries = [](Eigen::MatrixXd& d2DdXdP, const double& du,
                             const double& dv, const Vector3d& du_p,
                             const Vector3d& dv_p, const int& index_x,
                             const int& index_p) {
            d2DdXdP.block(index_x, index_p, 1, 3) =
                ((du_p * dv - du * dv_p) / (dv * dv)).transpose();
        };

        d2DdXdP.setZero();

        double du = (P2 - P1).transpose() * v;
        double dv = v_norm;
        Vector3d du_p, dv_p;
        //// G[0] = (P2 - P1).transpose() * (v / v_norm);
        {  //--- pG0pP1
            du_p = -1.0 * I * v + (1.0 - X[0]) * I * (P2 - P1);
            dv_p = (1.0 - X[0]) * v / v_norm;
            setEntries(d2DdXdP, du, dv, du_p, dv_p, 0, 0);
        }
        {  //--- pG0pP2
            du_p = I * v + X[0] * I * (P2 - P1);
            dv_p = X[0] * v / v_norm;
            setEntries(d2DdXdP, du, dv, du_p, dv_p, 0, 3);
        }

        //// G[1] = distanceWeight * (P4 - P3).transpose() * (v / v_norm);
        du = -1.0 * (P4 - P3).transpose() * v;
        {  //--- pG1pP1
            du_p = -(1.0 - X[0]) * I * (P4 - P3);
            dv_p = (1.0 - X[0]) * v / v_norm;
            setEntries(d2DdXdP, du, dv, du_p, dv_p, 1, 0);
        }
        {  //--- pG1pP2
            du_p = -X[0] * I * (P4 - P3);
            dv_p = X[0] * v / v_norm;
            setEntries(d2DdXdP, du, dv, du_p, dv_p, 1, 3);
        }
    }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
double Capsule::compute_D(const primitive_t& other) const { return -1; }

void Capsule::compute_dDdP(VectorXd& dDdP, const primitive_t& other) const {}

void Capsule::compute_d2DdP2(MatrixXd& d2DdP2, const primitive_t& other) const {
}
void Capsule::compute_d2DdP2_other(MatrixXd& d2DdP2_other,
                                   const primitive_t& other) const {}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
}  // namespace DCA

#endif /* __DCA_CAPSULE_H__ */