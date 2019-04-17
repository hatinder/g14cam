//
// Created by hsingh9 on 12/04/2019.
//

#include <fstream>
#include <iomanip>
#include <iostream>
#include "GaussQuadrature.hpp"

map<double, double> GaussQuadrature::LegendrePolynomial (int n, vector<double> v)
{
    int k = n - 1;
    map<double, double> results;    //C++ Dictionary map with key = x and value = calc results
    if (n == 0)
    {
        for (double i : v)
        {
            results[i] = 1;    //n=0 f(x)=1
        }
    }
    else if (n == 1)
    {
        for (double i : v)
        {
            results[i] = i; //n=1 f(x)=x
        }
    }
    else if (n > 1)
    {
        map<double, double> PhiK = LegendrePolynomial(k, v);        //recursive call to get Phi(k)(x)
        map<double, double> PhiK_1 = LegendrePolynomial(k - 1, v);  //recursive call to get Phi(k-1)(x)
        for (double i : v)                          //iterate thru all values of x
        {
            double firstTerm = (2 * k + 1) * i * PhiK[i];
            double secondTerm = (k) * PhiK_1[i];
            double mFactor = 1.0 / (k + 1);
            results[i] = mFactor * (firstTerm - secondTerm); //results to map with key=x and value=results
        }
    }
    return results;
}

map<double, double> GaussQuadrature::LegendrePolynomialDerivative (int n, vector<double> v)
{
    int k = n - 1;
    map<double, double> results;
    if (n == 0)
    {
        for (double i : v)
        {
            results[i] = 0;    //n=0 f'(x)=0
        }
    }
    else if (n == 1)
    {
        for (double i : v)
        {
            results[i] = 1; //n=1 f'(x)=1
        }
    }
    else if (n > 1)
    {
        map<double, double> PhiK = LegendrePolynomial(k, v);                    //recursive call to get Phi(k)(x)
        map<double, double> dPhiK = LegendrePolynomialDerivative(k, v);         //recursive call to get Phi'(k)(x)
        map<double, double> dPhiK_1 = LegendrePolynomialDerivative(k - 1, v);   //recursive call to get Phi'(k-1)(x)
        for (double i : v)
        {
            double firstTerm = (2 * k + 1) * (i * dPhiK[i] + PhiK[i]);
            double secondTerm = (k) * dPhiK_1[i];
            double mFactor = 1.0 / (k + 1);
            results[i] = mFactor * (firstTerm - secondTerm); //results to map with key=x and value=results
        }
    }
    return results;
}

//void
//GaussQuadrature::writeToFile (const string fNamePrefix, map<double, double> uniEvalPoints, const int k,
//                              vector<string> colNames)
//{
//
//}