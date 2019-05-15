//
// Created by hsingh9 on 22/04/2019.
//

#include <cmath>
#include <iostream>
#include "ODE.hpp"
#include "helper.hpp"

vector<map<double, double>>
ODE::applyRungeKutta2 (double (*f1) (double, double, double), double (*f2) (double), double u0, double v0, double T,
                       double N)
{
    vector<map<double, double>> keyTnValueYnPair;    //t_n and y_n key value pairs as return type for u and v
    double dt = T / N;
    double t0 = 0.0;
    map<double, double> vKeyValuePair;
    map<double, double> uKeyValuePair;
    vKeyValuePair[t0] = v0;
    uKeyValuePair[t0] = u0;
    for (int i = 0; i < N; ++i)
    {
        double vk1 = dt * (f1(v0, u0, t0));
        double uk1 = dt * (f2(v0));
        double vk2 = dt * (f1(v0 + vk1, u0 + uk1, t0 + dt));
//        double vk2 = dt * (f1(v0 + vk1, u0, t0));
//        double vNext = v0 + (vk1 + vk2) / 2.0;              //vNext = v(t+dt)
        double uk2 = dt * (f2(v0 + vk1));
        u0 += (uk1 + uk2) / 2.0;
//        v0 = vNext;
        v0 += (vk1 + vk2) / 2.0;
        t0 += dt;
        vKeyValuePair[t0] = v0;
        uKeyValuePair[t0] = u0;
    }
    keyTnValueYnPair.push_back(vKeyValuePair);
    keyTnValueYnPair.push_back(uKeyValuePair);
    return keyTnValueYnPair;
}

vector<map<double, double>>
ODE::applyImplicitMidpoint (double (*f1) (double, double, double, double, double), double (*f2) (double, double),
                            double u0, double v0, double T, double N)
{
    vector<map<double, double>> keyTnValueYnPair;    //t_n and y_n key value pairs as return type for u and v
    double dt = T / N;
    double t0 = 0.0;
    map<double, double> vKeyValuePair;
    map<double, double> uKeyValuePair;
    vKeyValuePair[t0] = v0;
    uKeyValuePair[t0] = u0;
    for (int i = 0; i < N; ++i)
    {
        double vk1Initial = v0 + dt * (f1(v0, u0, t0, v0, u0)), vk1 = 0.0;
        double uk1Initial = u0 + dt * (f2(v0, vk1Initial)), uk1 = 0.0;
        double TOL = pow(10, -12), tolerance = 1.0;
        int j = 0;
        for (; j < 100 && tolerance > TOL; ++j)
        {
            vk1 = v0 + dt * (f1(v0, u0, t0 + dt / 2.0, vk1Initial, uk1Initial));
            uk1 = u0 + dt * (f2(v0, vk1));
            tolerance = abs(uk1 - uk1Initial);
            vk1Initial = vk1;
            uk1Initial = uk1;
        }
//        cout<<"Iteration: "<<j<<endl;
        v0 = vk1;
        u0 = uk1;
        t0 += dt;
        vKeyValuePair[t0] = v0;
        uKeyValuePair[t0] = u0;
    }
    keyTnValueYnPair.push_back(vKeyValuePair);
    keyTnValueYnPair.push_back(uKeyValuePair);
    return keyTnValueYnPair;
}

vector<map<double, double>>
ODE::applyRungeKutta2 (vector<double, allocator<double>> (*F) (double, double, double, double), vector<double> initVal,
                       double T, double N)
{
    vector<map<double, double>> keyValuePair;
    double dt = T / N, t0 = 0.0;
    map<double, double> theta1Pair, theta2Pair, p1Pair, p2Pair;
    theta1Pair[t0] = initVal[0];
    theta2Pair[t0] = initVal[1];
    p1Pair[t0] = initVal[2];
    p2Pair[t0] = initVal[3];
    for (int i = 0; i < N; ++i)
    {
        vector<double> K1 = dt * (F(initVal[0], initVal[1], initVal[2],
                                    initVal[3]));  //helper.hpp has overloaded scalar product for vector
        vector<double> K2 = dt * (F(initVal[0] + K1[0], initVal[1] + K1[1], initVal[2] + K1[2], initVal[3] + K1[3]));
        double t = 2.0;
        initVal = initVal + (K1 + K2) / t;    // overloaded operators in helper.hpp
        t0 += dt;
        theta1Pair[t0] = initVal[0];
        theta2Pair[t0] = initVal[1];
        p1Pair[t0] = initVal[2];
        p2Pair[t0] = initVal[3];
    }
    keyValuePair.push_back(theta1Pair);
    keyValuePair.push_back(theta2Pair);
    keyValuePair.push_back(p1Pair);
    keyValuePair.push_back(p2Pair);
    return keyValuePair;
}

vector<map<double, double>>
ODE::applyImplicitMidpoint (vector<double, allocator<double>> (*F) (double, double, double, double),
                            vector<double> initVal, double T, double N)
{
    vector<double> Y(initVal.size(), 0.0);
    vector<map<double, double>> keyValuePair;
    double dt = T / N, t0 = 0.0, two = 2.0;
    map<double, double> theta1Pair, theta2Pair, p1Pair, p2Pair;
    theta1Pair[t0] = initVal[0];
    theta2Pair[t0] = initVal[1];
    p1Pair[t0] = initVal[2];
    p2Pair[t0] = initVal[3];
    for (int i = 0; i < N; ++i)
    {
        vector<double> yInitial = initVal + dt * F(initVal[0], initVal[1], initVal[2], initVal[3]);
        double TOL = pow(10, -12), tolerance = 1.0;
        int j = 0;
        for (; j < 100 && tolerance > TOL; ++j)
        {
            vector<double> mpValue = (yInitial + initVal) / two;
            Y = initVal + dt * F(mpValue[0], mpValue[1], mpValue[2], mpValue[3]);
            tolerance = abs(l2_norm(Y - yInitial));
            yInitial = Y;
        }
//        cout<<"Iteration: "<<j<<endl;
        initVal = Y;
        t0 += dt;
        theta1Pair[t0] = initVal[0];
        theta2Pair[t0] = initVal[1];
        p1Pair[t0] = initVal[2];
        p2Pair[t0] = initVal[3];
    }
    keyValuePair.push_back(theta1Pair);
    keyValuePair.push_back(theta2Pair);
    keyValuePair.push_back(p1Pair);
    keyValuePair.push_back(p2Pair);
    return keyValuePair;
}

map<double,vector<double>>
ODE::applyRK2(vector<double> (*F)(double, vector<double>), vector<double> initVal, double T, double N) {
    map<double,vector<double>> keyTnValueYnPair;    //t_n and y_n key value pairs as return type for u and v
    double dt = T / N;
    double t0 = 0.0;
    vector<double> yVal=initVal;
    keyTnValueYnPair[t0] = yVal;
    for (int i = 0; i < N; ++i)
    {
        double Dr=2.0;  // double data type for operator overload
        vector<double> k1 = dt * F(t0,yVal);
        vector<double> k2 = dt * F(t0+dt,yVal+k1);
        yVal=yVal + (k1+k2)/Dr;
        t0 += dt;
        keyTnValueYnPair[t0] = yVal;
    }
    return keyTnValueYnPair;
}

/*
vector<map<double, double>>
ODE::applyRungeKutta2 (vector<double, allocator<double>> (*F) (double, double, double, double), vector<double> initVal,
                       double T, double N)
{
    vector<map<double, double>> keyValuePair;
    double dt = T / N, t0=0.0;
    map<double, double> theta1Pair,theta2Pair,p1Pair,p2Pair;
    theta1Pair[t0] = initVal[0];
    theta2Pair[t0] = initVal[1];
    p1Pair[t0] = initVal[2];
    p2Pair[t0] = initVal[3];
    for (int i = 0; i < N; ++i)
    {
        vector<double> K1 = dt * (F(initVal[0], initVal[1], initVal[2],
                                    initVal[3]));  //helper.hpp has overloaded scalar product for vector
        vector<double> K2 = dt * (F(initVal[0] + K1[0], initVal[1] + K1[1], initVal[2] + K1[2], initVal[3] + K1[3]));
        double t=2.0;
        initVal = initVal + (K1+K2)/t;    // overloaded operators in helper.hpp
        t0 +=dt;
        theta1Pair[t0] = initVal[0];
        theta2Pair[t0] = initVal[1];
        p1Pair[t0] = initVal[2];
        p2Pair[t0] = initVal[3];
    }
    keyValuePair.push_back(theta1Pair);
    keyValuePair.push_back(theta2Pair);
    keyValuePair.push_back(p1Pair);
    keyValuePair.push_back(p2Pair);
    return keyValuePair;
}
*/


map<double,vector<double>>
ODE::applyIMP (vector<double> (*F) (double, vector<double>), vector<double> initVal, double T, double N)
{
    map<double,vector<double>> keyTnValueYnPair;    //t_n and y_n key value pairs as return type for u and v
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
