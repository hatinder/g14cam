#include <iostream>
#include <Dense>
#include <cmath>
#include "CubicSpline.hpp"

using namespace std;
using namespace Eigen;
int main ()
{
    std::cout << "Hello, World!" << std::endl;

    const int n=128;
    const double a=0,b=1;
    VectorXd nodalPoints=VectorXd::LinSpaced(n+1,a,b);
    VectorXd fValue(n+1);
    for (int i = 0; i < n+1; ++i)
    {
        fValue[i]=cos(6*M_PI*nodalPoints[i]);
    }
    cout<<"Function Value"<<endl;
    cout<<fValue<<endl;
    cout<<"Function Value"<<endl;

    CubicSpline cs;
    VectorXd coefficents=cs.findCoefficients(fValue);
    cout<<coefficents<<endl;
    cs.findSplineValues(coefficents,nodalPoints);
    return 0;
}