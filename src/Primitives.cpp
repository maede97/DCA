#include <DCA/Capsule.h>
#include <DCA/Primitives.h>
#include <DCA/Sphere.h>

namespace DCA {
double compute_D(const pair_t &pair,
                 const std::vector<primitive_t> &primitives) {
    return std::visit([](const auto &cp1,
                         const primitive_t &cp2) { return cp1.compute_D(cp2); },
                      primitives.at(pair.first), primitives.at(pair.second));
}

/**
 * Computes the gradient of a given pair and primitives
 */
void compute_dDdP(Eigen::VectorXd &dDdP, const pair_t &pair,
                  const std::vector<primitive_t> &primitives) {
    std::visit([&dDdP](const auto &cp1,
                       const primitive_t &cp2) { cp1.compute_dDdP(dDdP, cp2); },
               primitives.at(pair.first), primitives.at(pair.second));
}

/**
 * Computes the hessian of a given pair and primitives
 */
void compute_d2DdP2(Eigen::MatrixXd &d2DdP2, const pair_t &pair,
                    const std::vector<primitive_t> &primitives) {
    std::visit(
        [&d2DdP2](const auto &cp1, const primitive_t &cp2) {
            cp1.compute_d2DdP2(d2DdP2, cp2);
        },
        primitives.at(pair.first), primitives.at(pair.second));
}

void getParameter(VectorXd &parameter, const primitive_t primitive) {
    std::visit(
        overloaded{
            [&parameter](const Sphere &s) { s.getParameter(parameter); },
            [&parameter](const Capsule &c) { c.getParameter(parameter); }},
        primitive);
}

}  // namespace DCA