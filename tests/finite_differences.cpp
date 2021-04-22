#include <DCA/DCA.h>

#include <iostream>

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

void getParameters(Eigen::VectorXd& params, const DCA::primitive_t& prim) {
    std::visit(
        DCA::overloaded{
            [&params](const DCA::Sphere& s) { s.getParameters(params); },
            [&params](const DCA::Capsule& c) { c.getParameters(params); }},
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
        getParameters(params2, p2);

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
        getParameters(params2, p2);

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

int main(int argc, char const* argv[]) {
    /**
     * Run FD checks for each pair.
     */

    DCA::Sphere p1(Vector3d(1, 4., 2.), 0.36);
    DCA::Sphere p2(Vector3d(3., 6., 1.), 0.4);
    DCA::Capsule p3(Vector3d(0, 0, 0), Vector3d(1., 0, 0), 0.64);
    DCA::Capsule p4(Vector3d(2.4, 3., 0), Vector3d(3., 5., 1.), 0.5);

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

    return 0;
}