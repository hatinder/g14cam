//
// Created by hsingh9 on 12/04/2019.
//

#ifndef CODE_GAUSSQUADRATURE_HPP
#define CODE_GAUSSQUADRATURE_HPP


#include "IGaussQuadrature.hpp"
#include <string>
#include <vector>

using namespace std;

class GaussQuadrature : public IGaussQuadrature
{
public:

    map<double, double> LegendrePolynomial (int n, vector<double> v) override;

    map<double, double> LegendrePolynomialDerivative (int n, vector<double> v) override;

};


#endif //CODE_GAUSSQUADRATURE_HPP
