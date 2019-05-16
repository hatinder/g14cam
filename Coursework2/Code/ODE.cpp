//
// Created by hsingh9 on 22/04/2019.
//

#include <cmath>
#include <iostream>
#include "ODE.hpp"
#include "helper.hpp"


/*  Runge Kutta 2
    Algorithm is inline with Problem 2 (Explicit) Runge-Kutta 2 Method.
    Operator overloaded functions are in helper.cpp
    2a - It returns time as key value and x1 & x2 as u and u'
    2b,2c - It returns time as key value and theta1, theta2, p1 & p2
 */
map<double,vector<double>>
ODE::applyRK2(vector<double> (*F)(double, vector<double>), vector<double> initVal, double T, double N) {
    map<double,vector<double>> keyTnValueYnPair;
    double dt = T / N;
    double t0 = 0.0;
    vector<double> yVal=initVal;
    keyTnValueYnPair[t0] = yVal;
    for (int i = 0; i < N; ++i)
    {
        double Dr=2.0;  // double data type for operator overload
        vector<double> k1 = dt * F(t0,yVal);    // * operator overload
        vector<double> k2 = dt * F(t0+dt,yVal+k1); // * + operator overload
        yVal=yVal + (k1+k2)/Dr; // + / operator overload
        t0 += dt;
        keyTnValueYnPair[t0] = yVal;
    }
    return keyTnValueYnPair;
}

/*   Implicit Midpoint
    Algorithm is inline with Problem 2 Implicit Midpoint Rule.
    Operator overloaded functions are in helper.cpp
    2a - It returns time as key value and x1 & x2 as u and u'
    2b,2c - It returns time as key value and theta1, theta2, p1 & p2
 */
map<double,vector<double>>
ODE::applyIMP (vector<double> (*F) (double, vector<double>), vector<double> initVal, double T, double N)
{
    map<double,vector<double>> keyTnValueYnPair;
    double dt = T / N;
    double t0 = 0.0;
    double Dr=2.0;  // double data type for operator overload
    vector<double> yVal=initVal;
    keyTnValueYnPair[t0] = yVal;
    for (int i = 0; i < N; ++i)
    {
        yVal=initVal + dt*F(t0+dt/Dr,(yVal+initVal)/Dr);
        double TOL = pow(10, -12), tolerance = 1.0;
        int j = 0;
        for (; j < 100 && tolerance > TOL; ++j)
        {
            vector<double> yNext=initVal + dt*F(t0+dt/Dr,(yVal+initVal)/Dr);
            tolerance=l2_norm(yNext-yVal);
            yVal=yNext;
        }
        if(j==100)
        {
            cout<<"Warning: Problem Converging!!"<<endl;
            cout<<"Iteration: "<<j<<" , T:"<<T<<" , N: "<<N<<endl;
        }
        t0 += dt;
        keyTnValueYnPair[t0] = yVal;
        initVal=yVal;
    }
    return keyTnValueYnPair;
}
