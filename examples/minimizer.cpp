/**
 * This file should serve as an example on how to use the newton minimizer.
 * 
 * @author: Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#include "DCA/dca.h"

using namespace DCA;

int main(int argc, char const* argv[]) {
    CapsuleDistanceObjective f;

    NewtonOptimizer minimizer;

    Eigen::VectorXd x;
    x.resize(2);
    x << 0.5, 0.5;

    Eigen::VectorXd P_vector;
    P_vector.resize(12);
    P_vector << 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 2;

    minimizer.optimize(f, P_vector, x, 100);
    std::cout << x << std::endl;

    std::cout << "Distance: " << f.compute_O(P_vector, x) << std::endl;

    return 0;
}
