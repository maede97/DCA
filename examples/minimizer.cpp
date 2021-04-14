/**
 * This file should serve as an example on how to use the newton minimizer.
 * 
 * @author: Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#include "DCA/DCA.h"

using namespace DCA;

int main(int argc, char const* argv[]) {
    CapsuleDistanceObjective f;

    NewtonOptimizer minimizer;

    Eigen::VectorXd X;
    X.resize(2);
    X << 0.5, 0.5;

    Eigen::VectorXd P_vector;
    P_vector.resize(12);
    P_vector << Vector3d(0,0,0), Vector3d(1, 0, 0), Vector3d(0, 1, 0), Vector3d(1, 2, 0);

    bool converged =  minimizer.optimize(f, P_vector, X, 1000);
    std::cout << "Solver converged: " << (converged ? "true" : "false") << std::endl;
    std::cout << "Resulting X: " << X.transpose() << std::endl;
    std::cout << "Distance: " << f.compute_O(P_vector, X) << std::endl;

    return 0;
}
