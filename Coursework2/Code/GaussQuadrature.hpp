//
// Created by hsingh9 on 12/04/2019.
//

#ifndef CODE_GAUSSQUADRATURE_HPP
#define CODE_GAUSSQUADRATURE_HPP


#include "IGaussQuadrature.hpp"
#include <string>
#include <vector>
#include <cmath>


using namespace std;

class GaussQuadrature : public IGaussQuadrature
{
public:
    double functionValue(double x, double y){return sin(x*x+y+y);}
    map<double, double> LegendrePolynomial (int n, vector<double> v) override;

    map<double, double> LegendrePolynomialDerivative (int n, vector<double> v) override;

    vector<double> NewtonMethod (int n) override;

    vector<vector<double>> findPointsAndWeights (int n) override;

    double getExactValueFor1D (double a, double b) override;

    inline double getExactValueFor1C (double a, double b, int n) override {return ((pow(b,n+1)/(n+1))-(pow(a,n+1)/(n+1)));};

    double getApproxValue1C (double a, double b, vector<double> roots, vector<double> weights, int n, int degree) override;

    double
    getApproxValue1D (double a, double b, vector<double> roots, vector<double> weights, int n) override;

    double getApproxValue1E (double a, double b, double (*f) (double, double), int n, vector<double> roots,
                             vector<double> weights) override;

    double getApproxValue1F (double a, double b, double (*f) (double, double), int n, vector<double> roots,
                             vector<double> weights) override;
};


#endif //CODE_GAUSSQUADRATURE_HPP
