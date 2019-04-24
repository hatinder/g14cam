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
    virtual map<double, double>
    LegendrePolynomial (int n, vector<double> v) = 0;     // Legendre Polynomial Virtual Function
    virtual map<double, double>
    LegendrePolynomialDerivative (int n, vector<double> v) = 0;     // Legendre Polynomial Virtual Function
    virtual vector<double> NewtonMethod (int n) = 0;

    virtual vector<vector<double>> findPointsAndWeights (int n) = 0;

    virtual double getExactValueFor1C (double a, double b, int n) = 0;

    virtual double getExactValueFor1D (double a, double b) = 0;

    virtual double
    getApproxValue1C (double a, double b, vector<double> roots, vector<double> weights, int n, int degree) = 0;

    virtual double
    getApproxValue1D (double a, double b, vector<double> roots, vector<double> weights, int n) = 0;

    virtual double
    getApproxValue1E (double a, double b, double (*f) (double, double), int n, vector<double> roots,
                      vector<double> weights) = 0;

    virtual double
    getApproxValue1F (double a, double b, double (*f) (double, double), int n, vector<double> roots,
                      vector<double> weights) = 0;

};

#endif //CODE_IGAUSSQUADRATURE_HPP
