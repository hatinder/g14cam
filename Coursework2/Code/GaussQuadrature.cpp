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

vector<double> GaussQuadrature::NewtonMethod (int n)
{
    int maxIter = 100, m, k = 0;
    double TOL = pow(10, -12);
    vector<double> results(n, 0.0);
    if (n % 2 == 1)      // Exploit the fact: for n = odd number one root is zero. n=k+1
    {
        results[n - 1] = 0;
        m = n - 1;
    }
    else
    {
        m = n;
    }
    for (int i = 0; i < m / 2; i++)
    {
        vector<double> initialGuess(1, 0.0), nextValue(1, 0.0);
        double ig = cos(M_PI * (i + 1.0 - 1.0 / 4.0) / (n + 1.0 / 2.0));
        initialGuess[0] = ig;
        int j = 0;
        double tolerance = 1.0;
        for (; j < maxIter && tolerance > TOL; ++j)
        {
            auto Nr = LegendrePolynomial(n, initialGuess);             //Numerator = Phi(x)
            auto Dr = LegendrePolynomialDerivative(n, initialGuess);   //Denominator = dPhi(x) derivative
            nextValue[0] = initialGuess[0] - (Nr[initialGuess[0]] / Dr[initialGuess[0]]);
            tolerance = abs(nextValue[0] - initialGuess[0]);
            initialGuess[0] = nextValue[0];
        }
//        cout<<"Iteration: "<<j<<endl;
        results[k] = initialGuess[0];
        results[k + 1] = -initialGuess[0];  //Exploit the fact: roots appear symmetrically about x=0
        k = k + 2;
    }
    return results;
}

vector<vector<double>> GaussQuadrature::findPointsAndWeights (int n)
{
    vector<vector<double>> pointsNWeights;
    vector<double> roots(n, 0.0), weights(n, 0.0);     //intialize roots and weights with value 0.0 and size = n
    roots = NewtonMethod(n);                          //get all roots first
    map<double, double> dPhiK = LegendrePolynomialDerivative(n, roots); //get all derivatives as key value pairs
    for (int i = 0; i < n; ++i)                     //iterate thru roots to get corresponding weights
    {
        weights[i] = 2.0 / ((1 - roots[i] * roots[i]) * (dPhiK[roots[i]] * dPhiK[roots[i]]));
    }
    pointsNWeights.push_back(roots);    //first vector = roots
    pointsNWeights.push_back(weights);  //second vector = weights
    return pointsNWeights;
}

double
GaussQuadrature::getApproxValue1C (double a, double b, vector<double> roots, vector<double> weights, int n, int degree) //generic solution
{
    double result=0.0;
    for (int i = 0; i < n; ++i)
    {
        result+=pow((((b-a)*roots[i]+(b+a))/2.0),degree)*weights[i];
    }
    result=((b-a)/2.0)*result;
    return result;
}

double GaussQuadrature::getExactValueFor1D (double a, double b)
{
    return ((pow(sin(b),3)/3.0)-(pow(sin(a),3)/3.0));
}

double
GaussQuadrature::getApproxValue1D (double a, double b, vector<double> roots, vector<double> weights, int n)
{
    double result=0.0;
    for (int i = 0; i < n; ++i)
    {
        result+=sin(((b-a)*roots[i]+(b+a))/2.0)*sin(((b-a)*roots[i]+(b+a))/2.0)*cos(((b-a)*roots[i]+(b+a))/2.0)*weights[i];
    }
    result=((b-a)/2.0)*result;
    return result;
}

double GaussQuadrature::getApproxValue1E (double a, double b, double (*f) (double, double), int n, vector<double> roots,
                                          vector<double> weights)
{
    double x, y;
    double result=0.0;
    for (int i = 0; i < n; ++i)
    {
        y=((b-a)*roots[i]+b+a)/2.0;
        for (int j = 0; j < n; ++j)
        {
            x=((b-a)*roots[j]+b+a)/2.0;
            double fValue=f(x,y);
            result+=weights[j]*fValue*weights[i];
        }
    }
    return ((b-a)/2.0)*((b-a)/2.0)*result;
}

double GaussQuadrature::getApproxValue1F (double a, double b, double (*f) (double, double), int n, vector<double> roots,
                                          vector<double> weights)
{
    double x, y;
    double result=0.0;
    for (int i = 0; i < n; ++i)
    {
        x=((b-a)*roots[i]+b+a)/2.0;
        double yout=0.0;
        for (int j = 0; j < n; ++j)
        {
            y=((x-a)*roots[j]+x+a)/2.0;
            double fValue=f(x,y);
            yout+=weights[j]*fValue;
        }
        result+=((x-a)/2.0)*yout*weights[i];
    }
    return ((b-a)/2.0)*result;
}

