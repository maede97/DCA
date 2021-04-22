#ifndef __DCA_AD_POINTTOLINESEGMENT_DISTANCE_H__
#define __DCA_AD_POINTTOLINESEGMENT_DISTANCE_H__

#include <Eigen/Core>

namespace PointToLineSegmentDistance_CodeGen {
void AD_PointToLineSegmentDistance(const Eigen::Matrix<double, 9, 1>& x_e,
                                   double sigScale, double& objVal);

void AD_PointToLineSegmentDistanceGradient(
    const Eigen::Matrix<double, 9, 1>& x_e, double sigScale,
    Eigen::Matrix<double, 9, 1>& gradient);

void AD_PointToLineSegmentDistanceHessian(
    const Eigen::Matrix<double, 9, 1>& x_e, double sigScale,
    Eigen::Matrix<double, 9, 9>& hessian);

}  // namespace PointToLineSegmentDistance_CodeGen

#endif /* __DCA_AD_POINTTOLINESEGMENT_DISTANCE_H__ */