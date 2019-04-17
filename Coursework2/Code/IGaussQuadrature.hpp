//
// Created by hsingh9 on 16/04/2019.
//

#ifndef CODE_IGAUSSQUADRATURE_HPP
#define CODE_IGAUSSQUADRATURE_HPP

#include <vector>
#include <map>

using namespace std;


class IGaussQuadrature
{
public:
    virtual map<double, double> LegendrePolynomial (int n, vector<double> v) = 0;     // Legendre Polynomial Virtual Function
    virtual map<double, double> LegendrePolynomialDerivative (int n, vector<double> v) = 0;     // Legendre Polynomial Virtual Function
};

#endif //CODE_IGAUSSQUADRATURE_HPP
