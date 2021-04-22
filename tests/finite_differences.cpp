#include <DCA/DCA.h>

#include <iostream>

// #define RUN_INDIVIDUAL_CAPSULE_TEST

using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::VectorXd;

unsigned int getParameterSize(const DCA::primitive_t& prim) {
    return std::visit(
        DCA::overloaded{
            [](const DCA::Sphere& s) { return s.getParametersSize(); },
            [](const DCA::Capsule& c) { return c.getParametersSize(); }},
        prim);
}

const char* getPrimitiveType(DCA::primitive_t prim) {
    return std::visit(
        DCA::overloaded{[](const DCA::Sphere& s) { return "Sphere"; },
                        [](const DCA::Capsule& s) { return "Capsule"; }},
        prim);
}

DCA::primitive_t makeNewFromParameters(const DCA::primitive_t& old,
                                       const Eigen::VectorXd& params) {
    const char* name = getPrimitiveType(old);

    if (name == "Sphere") {
        return DCA::Sphere(params);
    } else if (name == "Capsule") {
        return DCA::Capsule(params);
    }
    LOG << ANSI_COLOR_RED << "Could not create a new primitive from old.";
    throw std::logic_error("");
}

void check_dDdP_FD_entries(const VectorXd& dDdP, const VectorXd& dDdP_FD) {
    if (dDdP.size() != dDdP_FD.size()) {
        LOG << ANSI_COLOR_RED << "SIZE DOES NOT MATCH!";
        return;
    }

    std::cout << "Analytical norm: " << dDdP.norm()
              << ", FD norm: " << dDdP_FD.norm() << std::endl;

    bool hasMismatch = false;
    for (int i = 0; i < dDdP.size(); i++) {
        if (fabs(dDdP(i) - dDdP_FD(i)) > 1e-3) {
            std::cout << ANSI_COLOR_RED << "Mismatch element " << i << ": "
                      << dDdP(i) << " (analytical) vs. " << dDdP_FD(i)
                      << " (FD)" << ANSI_COLOR_DEFAULT << std::endl;
            hasMismatch = true;
        }
    }

    if (!hasMismatch)
        std::cout << ANSI_COLOR_GREEN << "  All Good." << ANSI_COLOR_DEFAULT
                  << std::endl;
}

void check_d2DdP2_FD_entries(const MatrixXd& d2DdP2,
                             const MatrixXd& d2DdP2_FD) {
    // Todo: check rows and cols
    if (d2DdP2.rows() != d2DdP2_FD.rows() ||
        d2DdP2.cols() != d2DdP2_FD.cols()) {
        LOG << ANSI_COLOR_RED << "SIZE DOES NOT MATCH!";
        return;
    }

    std::cout << "Analytical norm: " << d2DdP2.norm()
              << ", FD norm: " << d2DdP2_FD.norm() << std::endl;

    bool hasMismatch = false;
    for (int i = 0; i < d2DdP2.rows(); i++) {
        for (int j = 0; j < d2DdP2.cols(); j++) {
            if (fabs(d2DdP2(i, j) - d2DdP2_FD(i, j)) > 1e-3) {
                std::cout << ANSI_COLOR_RED << "Mismatch element " << i << "/"
                          << j << ": " << d2DdP2(i, j) << " (analytical) vs. "
                          << d2DdP2_FD(i, j) << " (FD)" << ANSI_COLOR_DEFAULT
                          << std::endl;
                hasMismatch = true;
            }
        }
    }

    if (!hasMismatch)
        std::cout << ANSI_COLOR_GREEN << "  All Good." << ANSI_COLOR_DEFAULT
                  << std::endl;
}

void check_dDdP_FD(DCA::primitive_t p1, DCA::primitive_t p2) {
    std::vector<DCA::primitive_t> customPrimitiveVector{p1, p2};
    DCA::pair_t customPair = {0, 1};

    // Now compute analytical gradient
    VectorXd dDdP;
    DCA::compute_dDdP(dDdP, customPair, customPrimitiveVector);

    VectorXd params2;

    VectorXd dDdP_fd;
    dDdP_fd.resize(dDdP.size());

    const double change = 0.0001;
    for (int i = 0; i < getParameterSize(p2); i++) {
        getParameter(params2, p2);

        // alter the params a bit
        params2(i) += change;

        // Alter the pair
        customPrimitiveVector.at(1) = makeNewFromParameters(p2, params2);

        double d1 = DCA::compute_D(customPair, customPrimitiveVector);

        params2(i) -= 2. * change;
        customPrimitiveVector.at(1) = makeNewFromParameters(p2, params2);
        double d2 = DCA::compute_D(customPair, customPrimitiveVector);

        // Compute the FD dDdP entry
        dDdP_fd(i) = (d1 - d2) / (2. * change);
    }

    std::cout << getPrimitiveType(p1) << " vs. " << getPrimitiveType(p2)
              << std::endl;
    check_dDdP_FD_entries(dDdP, dDdP_fd);
}

void check_d2DdP2_FD(DCA::primitive_t p1, DCA::primitive_t p2) {
    std::vector<DCA::primitive_t> customPrimitiveVector{p1, p2};
    DCA::pair_t customPair = {0, 1};

    MatrixXd d2DdP2;
    DCA::compute_d2DdP2(d2DdP2, customPair, customPrimitiveVector);

    VectorXd params2;

    MatrixXd d2DdP2_fd;
    d2DdP2_fd.resize(d2DdP2.rows(), d2DdP2.cols());

    const double change = 0.0001;
    for (int i = 0; i < getParameterSize(p2); i++) {
        getParameter(params2, p2);

        // alter the params a bit
        params2(i) += change;

        // Alter the pair
        customPrimitiveVector.at(1) = makeNewFromParameters(p2, params2);
        VectorXd dDdP1;
        DCA::compute_dDdP(dDdP1, customPair, customPrimitiveVector);

        params2(i) -= 2. * change;
        customPrimitiveVector.at(1) = makeNewFromParameters(p2, params2);
        VectorXd dDdP2;
        DCA::compute_dDdP(dDdP2, customPair, customPrimitiveVector);

        // Compute the FD dDdP entry
        d2DdP2_fd.col(i) = (dDdP1 - dDdP2) / (2. * change);
    }

    std::cout << getPrimitiveType(p1) << " vs. " << getPrimitiveType(p2)
              << std::endl;
    check_d2DdP2_FD_entries(d2DdP2, d2DdP2_fd);
}

#ifdef RUN_INDIVIDUAL_CAPSULE_TEST
/**
 * This function requires the private members in DCA::CapsuleDistanceHelper to be actually public.
 */
void check_d2DdP2_individual_FD(const DCA::Capsule& p1,
                                const DCA::Capsule& p2) {
    // Fancy helpers to make P const
    const VectorXd P = [&]() {
        VectorXd P_;
        P_.resize(12);
        VectorXd P1, P2;
        p2.getParameter(P2);
        p1.getParameter(P1);
        P_ << P2, P1;
        return P_;
    }();

    // Fancy helpers to make X const
    const VectorXd X = [&]() {
        VectorXd X_;
        DCA::CapsuleDistanceHelper::solveForX(P, X_);
        return X_;
    }();

    const double change = 0.0001;
    VectorXd params2;

    MatrixXd dXdP;  // size 2x6
    {
        MatrixXd dXdP_FD;
        DCA::CapsuleDistanceHelper::compute_dXdP(dXdP, P, X);

        VectorXd P_copy = P;

        dXdP_FD.resize(dXdP.rows(), dXdP.cols());

        for (int i = 0; i < p2.getParametersSize(); i++) {
            p2.getParameter(params2);

            params2(i) += change;
            P_copy.head(6) = params2;
            VectorXd X1;
            DCA::CapsuleDistanceHelper::solveForX(P_copy, X1);

            params2(i) -= 2. * change;
            P_copy.head(6) = params2;
            VectorXd X2;
            DCA::CapsuleDistanceHelper::solveForX(P_copy, X2);

            dXdP_FD.col(i) = (X1 - X2) / (2. * change);
        }

        std::cout << "dXdP" << std::endl;
        check_d2DdP2_FD_entries(dXdP, dXdP_FD);
    }

    MatrixXd p2DpX2;  // size 2x2
    {
        DCA::CapsuleDistanceHelper::objective().compute_d2DdX2(p2DpX2, P, X);
        MatrixXd p2DpX2_FD;
        p2DpX2_FD.resize(p2DpX2.rows(), p2DpX2.cols());

        for (int i = 0; i < X.size(); i++) {
            VectorXd dDdX1, dDdX2;

            params2 = X;

            params2(i) += change;
            DCA::CapsuleDistanceHelper::objective().compute_dDdX(dDdX1, P,
                                                                 params2);

            params2(i) -= 2. * change;
            DCA::CapsuleDistanceHelper::objective().compute_dDdX(dDdX2, P,
                                                                 params2);

            p2DpX2_FD.col(i) = (dDdX1 - dDdX2) / (2. * change);
        }

        std::cout << "p2DpX2" << std::endl;
        check_d2DdP2_FD_entries(p2DpX2, p2DpX2_FD);
    }

    MatrixXd p2DpP2;  // size 6x6
    {
        DCA::CapsuleDistanceHelper::objective().compute_d2DdP2(p2DpP2, P, X);
        VectorXd P_copy = P;

        MatrixXd p2DpP2_FD;
        p2DpP2_FD.resize(p2DpP2.rows(), p2DpP2.cols());

        VectorXd X_copy;
        for (int i = 0; i < p2.getParametersSize(); i++) {
            VectorXd dDdP1, dDdP2;

            p2.getParameter(params2);

            params2(i) += change;
            P_copy.head(6) = params2;
            DCA::CapsuleDistanceHelper::solveForX(P_copy, X_copy);
            DCA::CapsuleDistanceHelper::objective().compute_dDdP(dDdP1, P_copy,
                                                                 X_copy);

            params2(i) -= 2. * change;
            P_copy.head(6) = params2;
            DCA::CapsuleDistanceHelper::solveForX(P_copy, X_copy);
            DCA::CapsuleDistanceHelper::objective().compute_dDdP(dDdP2, P_copy,
                                                                 X_copy);

            p2DpP2_FD.col(i) = (dDdP1 - dDdP2) / (2. * change);
        }
        std::cout << "p2DpP2" << std::endl;
        check_d2DdP2_FD_entries(p2DpP2, p2DpP2_FD);
    }

    MatrixXd p2DpXpP;  // size 2x6

    {
        DCA::CapsuleDistanceHelper::objective().compute_d2DdXdP(p2DpXpP, P, X);
        VectorXd P_copy = P;
        MatrixXd p2DpXpP_FD;
        p2DpXpP_FD.resize(p2DpXpP.rows(), p2DpXpP.cols());

        for (int i = 0; i < p2.getParametersSize(); i++) {
            VectorXd dDdX1, dDdX2;
            VectorXd X_copy;

            p2.getParameter(params2);

            params2(i) += change;
            P_copy.head(6) = params2;
            DCA::CapsuleDistanceHelper::solveForX(P_copy, X_copy);
            DCA::CapsuleDistanceHelper::objective().compute_dDdX(dDdX1, P_copy,
                                                                 X_copy);

            params2(i) -= 2. * change;
            P_copy.head(6) = params2;
            DCA::CapsuleDistanceHelper::solveForX(P_copy, X_copy);
            DCA::CapsuleDistanceHelper::objective().compute_dDdX(dDdX2, P_copy,
                                                                 X_copy);

            p2DpXpP_FD.col(i) = (dDdX1 - dDdX2) / (2. * change);
        }

        std::cout << "p2DpXpP" << std::endl;
        check_d2DdP2_FD_entries(p2DpXpP, p2DpXpP_FD);
    }

    // size 6x6
    MatrixXd d2DdP2;
    d2DdP2 = dXdP.transpose() * p2DpX2 * dXdP + dXdP.transpose() * p2DpXpP +
             p2DpXpP.transpose() * dXdP + p2DpP2;
}

#endif

int main(int argc, char const* argv[]) {
    /**
     * Run FD checks for each pair.
     */

    DCA::Sphere p1(Vector3d(1, 4., 2.), 0.36);
    DCA::Sphere p2(Vector3d(3., 6., 1.), 0.4);
    DCA::Capsule p3(Vector3d(0, 0, 0), Vector3d(1., 0, 0), 0.64);
    DCA::Capsule p4(Vector3d(0.2, 0, 2), Vector3d(1.5, 0, 3), 0.5);

    std::cout << ANSI_COLOR_CYAN << "Checking dDdP" << ANSI_COLOR_DEFAULT
              << std::endl;
    check_dDdP_FD(p1, p2);  // Sphere - Sphere
    check_dDdP_FD(p1, p3);  // Sphere - Capsule
    check_dDdP_FD(p4, p2);  // Capsule - Sphere
    check_dDdP_FD(p3, p4);  // Capsule - Capsule

    std::cout << ANSI_COLOR_CYAN << "Checking d2DdP2" << ANSI_COLOR_DEFAULT
              << std::endl;
    check_d2DdP2_FD(p1, p2);  // Sphere - Sphere
    check_d2DdP2_FD(p1, p3);  // Sphere - Capsule
    check_d2DdP2_FD(p4, p2);  // Capsule - Sphere
    check_d2DdP2_FD(p3, p4);  // Capsule - Capsule

#ifdef RUN_INDIVIDUAL_CAPSULE_TEST
    std::cout << ANSI_COLOR_CYAN
              << "Checking d2DdP2 individual for Capsule - Capsule"
              << ANSI_COLOR_DEFAULT << std::endl;
    check_d2DdP2_individual_FD(p3, p4);
#endif

    return 0;
}