#include "DCA/dca.h"

using namespace DCA;

class CapsuleObjective : public Objective {
public:
    CapsuleObjective(Vector3d P1, Vector3d P2, Vector3d P3, Vector3d P4)
        : m_P1(P1), m_P2(P2), m_P3(P3), m_P4(P4) {}

    double evaluate(const VectorXd& x) const override {
        Vector3d P12 = m_P1 + x[0] * (m_P2 - m_P1);
        Vector3d P34 = m_P3 + x[1] * (m_P4 - m_P3);
        return (P12 - P34).norm();
    }

    VectorXd gradient(const VectorXd& x) const override {
        Vector3d P12 = m_P1 + x[0] * (m_P2 - m_P1);
        Vector3d P34 = m_P3 + x[1] * (m_P4 - m_P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;
        VectorXd pDpP;
        pDpP.resize(12);
        pDpP.segment(0, 3) = (1.0 - x[0]) * (v / v_norm);
        pDpP.segment(3, 3) = x[0] * (v / v_norm);
        pDpP.segment(6, 3) = -(1.0 - x[1]) * (v / v_norm);
        pDpP.segment(9, 3) = -x[1] * (v / v_norm);

        return pDpP;
    }

    SparseMatrixd hessian(const VectorXd& x) const override {
        const double t12 = x[0];
        const double t34 = x[1];

        Vector3d P12 = m_P1 + t12 * (m_P2 - m_P1);
        Vector3d P34 = m_P3 + t34 * (m_P4 - m_P3);
        Vector3d v = P12 - P34;
        double v_norm = v.norm();
        if (v_norm < 1e-8) v_norm = 1e-8;

        const Matrix3d I_vn = v_norm * Matrix3d::Identity();
        const Matrix3d v_T_v = v * (v / v_norm).transpose();
        const double v_norm_sq = v_norm * v_norm;

        double du_p, dv_p;
        Matrix3d ppD;
        double du = 1.0 - t12;

        MatrixXd p2DpP2(12, 12);
        ////// pDpP1 = (1.0 - t12) * (v / v_norm)
        {  //--- dP1
            du_p = (1.0 - t12) * (1.0 - t12);
            dv_p = (1.0 - t12);
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(0, 0, 3, 3) = ppD;
        }
        {  //--- dP2
            du_p = t12 * (1.0 - t12);
            dv_p = t12;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(3, 0, 3, 3) = ppD;
            p2DpP2.block(0, 3, 3, 3) = ppD;
        }
        {  //--- dP3
            du_p = -(1.0 - t34) * (1.0 - t12);
            dv_p = -(1.0 - t34);
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(6, 0, 3, 3) = ppD;
            p2DpP2.block(0, 6, 3, 3) = ppD;
        }
        {  //--- dP4
            du_p = -t34 * (1.0 - t12);
            dv_p = -t34;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(9, 0, 3, 3) = ppD;
            p2DpP2.block(0, 9, 3, 3) = ppD;
        }
        ////// pDpP2 = t12 * (v / v_norm)
        du = t12;
        {  //--- dP2
            du_p = t12 * t12;
            dv_p = t12;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(3, 3, 3, 3) = ppD;
        }
        {  //--- dP3
            du_p = t12 * -(1.0 - t34);
            dv_p = -(1.0 - t34);
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(6, 3, 3, 3) = ppD;
            p2DpP2.block(3, 6, 3, 3) = ppD;
        }
        {  //--- dP4
            du_p = t12 * -t34;
            dv_p = -t34;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(9, 3, 3, 3) = ppD;
            p2DpP2.block(3, 9, 3, 3) = ppD;
        }
        ////// pDpP2 = -(1.0 - t34) * (v / v_norm)
        du = -(1.0 - t34);
        {  //--- dP3
            du_p = -(1.0 - t34) * -(1.0 - t34);
            dv_p = -(1.0 - t34);
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(6, 6, 3, 3) = ppD;
            p2DpP2.block(6, 6, 3, 3) = ppD;
        }
        {  //--- dP4
            du_p = -t34 * -(1.0 - t34);
            dv_p = -t34;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(9, 6, 3, 3) = ppD;
            p2DpP2.block(6, 9, 3, 3) = ppD;
        }
        ////// pDpP2 = -t34 * (v / v_norm)
        du = -t34;
        {  //--- dP4
            du_p = -t34 * -t34;
            dv_p = -t34;
            ppD = (du_p * I_vn - du * dv_p * v_T_v) / v_norm_sq;
            p2DpP2.block(9, 9, 3, 3) = ppD;
        }
        return p2DpP2.sparseView(); // for now
    }

private:
    Vector3d m_P1, m_P2, m_P3, m_P4;
};

int main(int argc, char const* argv[]) {
    CapsuleObjective* f = new CapsuleObjective(Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0),
                      Vector3d(1, 1, 0));
    
    NewtonOptimizer minimizer;

    Eigen::VectorXd x;
    x.resize(2);
    x << 0.5, 0.5;

    std::cout << minimizer.optimize(f, x, 100) << std::endl;

    std::cout << x << std::endl;

    delete f;

    return 0;
}
