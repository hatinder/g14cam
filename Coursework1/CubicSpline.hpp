//
// Created by hsingh9 on 08/03/2019.
//

#ifndef COURSEWORK1_CUBICSPLINE_HPP
#define COURSEWORK1_CUBICSPLINE_HPP
#include <iostream>
#include <Dense>

using namespace Eigen;
using namespace std;

class CubicSpline
{

public:
    VectorXd findCoefficients(VectorXd fValue);
    ArrayXXd findSplineValues (VectorXd coeff, VectorXd nPoints);
    void writeToFile(string name, ArrayXXd uniEvalPoints );
private:
    double getBx (double xVal);
    VectorXd computeBi(VectorXd nPoints, double x);

};


#endif //COURSEWORK1_CUBICSPLINE_HPP
