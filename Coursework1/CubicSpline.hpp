//
// Created by hsingh9 on 08/03/2019.
//

#ifndef COURSEWORK1_CUBICSPLINE_HPP
#define COURSEWORK1_CUBICSPLINE_HPP
#include <vector>
#include <Dense>
using namespace Eigen;

class CubicSpline
{

public:
    VectorXd findCoefficients(VectorXd fValue);
    VectorXd findSplineValues (VectorXd coeff, VectorXd nPoints);
private:
    double getBx (double xVal);
    VectorXd computeBi(VectorXd nPoints, double x);

};


#endif //COURSEWORK1_CUBICSPLINE_HPP
