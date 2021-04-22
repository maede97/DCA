#include <DCA/Helpers.h>
#include <DCA/AD_PointToLineSegmentDistance.h>

namespace DCA {

double CapsuleDistanceObjective::compute_O(const VectorXd& P,
                                           const VectorXd& X) const {
    double value = 0.0;

    //--- Shortest distance
    value += compute_D(P, X);

    //--- Regularizer
    for (int i = 0; i < maxRegularizerIndex; i++)
        value += regularizerWeight * 0.5 * (X[i] - 0.5) * (X[i] - 0.5);

    //--- Constraint
    value += constraintWeight * sulc.compute_F(-X[0]);
    value += constraintWeight * sulc.compute_F(X[0] - 1.0);
    value += constraintWeight * sulc.compute_F(-X[1]);
    value += constraintWeight * sulc.compute_F(X[1] - 1.0);

    return value;
}

void CapsuleDistanceObjective::compute_dOdX(VectorXd& dOdX, const VectorXd& P,
                                            const VectorXd& X) const {
    //--- Shortest distance
    compute_dDdX(dOdX, P, X);

    //--- Regularizer
    for (int i = 0; i < maxRegularizerIndex; i++)
        dOdX[i] += regularizerWeight * (X[i] - 0.5);

    //--- Constraint
    dOdX[0] -= constraintWeight * sulc.compute_dFdX(-X[0]);
    dOdX[0] += constraintWeight * sulc.compute_dFdX(X[0] - 1.0);
    dOdX[1] -= constraintWeight * sulc.compute_dFdX(-X[1]);
    dOdX[1] += constraintWeight * sulc.compute_dFdX(X[1] - 1.0);
}

void CapsuleDistanceObjective::compute_d2OdX2(MatrixXd& d2OdX2,
                                              const VectorXd& P,
                                              const VectorXd& X) const {
    //--- Shortest distance
    compute_d2DdX2(d2OdX2, P, X);
    d2OdX2(1, 0) = 0.;

    //--- Regularizer
    for (int i = 0; i < maxRegularizerIndex; i++)
        d2OdX2(i, i) += regularizerWeight;

    //--- Constraint
    d2OdX2(0, 0) += constraintWeight * sulc.compute_d2FdX2(-X[0]);
    d2OdX2(0, 0) += constraintWeight * sulc.compute_d2FdX2(X[0] - 1.0);
    d2OdX2(1, 1) += constraintWeight * sulc.compute_d2FdX2(-X[1]);
    d2OdX2(1, 1) += constraintWeight * sulc.compute_d2FdX2(X[1] - 1.0);
}

void CapsuleDistanceObjective::compute_dDdX(VectorXd& dOdX, const VectorXd& P,
                                            const VectorXd& X) const {
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

void CapsuleDistanceObjective::compute_d2DdX2(MatrixXd& d2OdX2,
                                              const VectorXd& P,
                                              const VectorXd& X) const {
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

double CapsuleDistanceObjective::compute_D(const VectorXd& P,
                                           const VectorXd& X) const {
    Vector3d P12 = P.head(3) + X[0] * (P.segment(3, 3) - P.head(3));
    Vector3d P34 = P.segment(6, 3) + X[1] * (P.tail(3) - P.segment(6, 3));
    return (P12 - P34).norm();
}

void CapsuleDistanceObjective::compute_dDdP(VectorXd& dDdP, const VectorXd& P,
                                            const VectorXd& X) const {
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

void CapsuleDistanceObjective::compute_d2DdP2(MatrixXd& d2DdP2,
                                              const VectorXd& P,
                                              const VectorXd& X) const {
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

void CapsuleDistanceObjective::compute_d2DdXdP(MatrixXd& d2DdXdP,
                                               const VectorXd& P,
                                               const VectorXd& X) const {
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

CapsuleDistanceObjective CapsuleDistanceHelper::objective =
    CapsuleDistanceObjective();

void CapsuleDistanceHelper::compute_dXdP(MatrixXd& dXdP, const VectorXd& P,
                                         const VectorXd& X) {
    dXdP.resize(2, 6);
    MatrixXd d2OdX2;
    objective.compute_d2OdX2(d2OdX2, P, X);

    MatrixXd dGdX = d2OdX2;
    dGdX(0, 1) = dGdX(1, 0);  // is this still needed?

    MatrixXd dGdP;
    objective.compute_d2DdXdP(dGdP, P, X);

    dXdP = -1. * dGdX.inverse() * dGdP;
}

double CapsuleDistanceHelper::compute_D(const VectorXd& P, const double& r12,
                                        const double& r34) {
    VectorXd X;
    solveForX(P, X);
    return objective.compute_D(P, X) - (r12 + r34);
}

/**
     * Computes the second derivative of the distance with respect to
     * the points of the first capsule.
     */
void CapsuleDistanceHelper::compute_dDdP(VectorXd& dDdP, const VectorXd& P) {
    VectorXd X;
    solveForX(P, X);

    MatrixXd dXdP;  // size 2x6
    compute_dXdP(dXdP, P, X);

    VectorXd pDpX;  // size 2
    objective.compute_dDdX(pDpX, P, X);

    VectorXd pDpP;  // size 6
    objective.compute_dDdP(pDpP, P, X);

    dDdP = dXdP.transpose() * pDpX + pDpP;
}

/**
     * Computes the first derivative of the distance with respect to
     * the points of the first capsule.
     */
void CapsuleDistanceHelper::compute_d2DdP2(MatrixXd& d2DdP2,
                                           const VectorXd& P) {
    VectorXd X;
    solveForX(P, X);

    MatrixXd dXdP;  // size 2x6
    compute_dXdP(dXdP, P, X);

    MatrixXd p2DpX2;  // size 2x2
    objective.compute_d2DdX2(p2DpX2, P, X);

    MatrixXd p2DpP2;  // size 6x6
    objective.compute_d2DdP2(p2DpP2, P, X);

    MatrixXd p2DpXpP;  // size 2x6
    objective.compute_d2DdXdP(p2DpXpP, P, X);

    // size 6x6
    d2DdP2 = dXdP.transpose() * p2DpX2 * dXdP + dXdP.transpose() * p2DpXpP +
             p2DpXpP.transpose() * dXdP + p2DpP2;
}

void CapsuleDistanceHelper::solveForX(const VectorXd& P, VectorXd& X) {
    NewtonOptimizer optimizer;
    X.resize(2);
    X << 0.5, 0.5;
    optimizer.optimize(objective, P, X, 100);
}

void CapsuleDistanceHelper::computeClosestPointOnLines(Vector3d& P12,
                                                       Vector3d& P34,
                                                       const VectorXd& P) {
    VectorXd X;
    solveForX(P, X);
    Vector3d P1 = P.segment(0, 3);
    Vector3d P2 = P.segment(3, 3);
    Vector3d P3 = P.segment(6, 3);
    Vector3d P4 = P.segment(9, 3);

    P12 = P1 + X[0] * (P2 - P1);
    P34 = P3 + X[1] * (P4 - P3);
}

double SphereCapsuleDistanceHelper::compute_D(const VectorXd& P,
                                              const double& rSphere,
                                              const double& rCapsule) {
    Vector3d P0 = P.head(3);
    return (P0 - computeClosestPointOnLine(P)).norm() - (rSphere + rCapsule);
}

void SphereCapsuleDistanceHelper::compute_dDdP(VectorXd& dDdP,
                                               const VectorXd& P) {
    dDdP.resize(9);

    Eigen::Matrix<double, 9, 1> grad;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceGradient(
        P, sigScale(), grad);
    dDdP = grad;
}

void SphereCapsuleDistanceHelper::compute_d2DdP2(MatrixXd& d2DdP2,
                                                 const VectorXd& P) {
    d2DdP2.resize(9, 9);
    Eigen::Matrix<double, 9, 9> hess;
    PointToLineSegmentDistance_CodeGen::AD_PointToLineSegmentDistanceHessian(
        P, sigScale(), hess);
    d2DdP2 = hess;
}

Vector3d SphereCapsuleDistanceHelper::computeClosestPointOnLine(
    const VectorXd& P) {
    Vector3d P1 = P.segment(3, 3);
    Vector3d P2 = P.tail(3);
    return P1 + (P2 - P1) * compute_t(P);
}

double SphereCapsuleDistanceHelper::compute_t(const VectorXd& P) {
    Vector3d P0 = P.head(3);
    Vector3d P1 = P.segment(3, 3);
    Vector3d P2 = P.tail(3);
    double t = -1.0 * (((P1 - P0).dot(P2 - P1)) / (P2 - P1).squaredNorm());
    return sigmoid(t, sigScale());
}

double SphereCapsuleDistanceHelper::sigScale() {
    double s = 5.0;
    return s;
}

}  // namespace DCA