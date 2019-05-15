//
// Created by hsingh9 on 22/04/2019.
//

#ifndef CODE_ODE_HPP
#define CODE_ODE_HPP

#include <vector>
#include <map>

using namespace std;

class ODE
{
public:
    vector<map<double, double>>
    applyRungeKutta2 (double (*f1) (double, double, double), double (*f2) (double), double u0, double v0, double T,
                      double N);

    map<double,vector<double>>
    applyRK2 (vector<double> (*F) (double, vector<double>), vector<double> initVal, double T, double N);

    map<double,vector<double>>
    applyIMP (vector<double> (*F) (double, vector<double>), vector<double> initVal, double T, double N);

    vector<map<double, double>>
    applyImplicitMidpoint (double (*f1) (double, double, double, double, double), double (*f2) (double, double),
                           double u0, double v0, double T, double N);

    vector<map<double, double>>
    applyRungeKutta2 (vector<double, allocator<double>> (*F) (double, double, double, double), vector<double> initVal,
                      double T, double N);

    vector<map<double, double>>
    applyImplicitMidpoint (vector<double, allocator<double>> (*F) (double, double, double, double),
                           vector<double> initVal, double T, double N);
};


#endif //CODE_ODE_HPP
