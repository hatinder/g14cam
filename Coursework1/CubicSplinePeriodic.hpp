//
// Created by hsingh9 on 15/03/2019.
//

#ifndef COURSEWORK1_CUBICSPLINEPERIODIC_HPP
#define COURSEWORK1_CUBICSPLINEPERIODIC_HPP

#include "CubicSpline.hpp"

class CubicSplinePeriodic : public CubicSpline
{
public:
    VectorXd findCoefficients(VectorXd fValue);
    ArrayXXd findSplineValues (VectorXd coeff, VectorXd nPoints);

protected:
    VectorXd computeBi(VectorXd nPoints, double x);
};


#endif //COURSEWORK1_CUBICSPLINEPERIODIC_HPP
